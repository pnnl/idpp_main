"""
    idpp/db/util.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with misc utilities for IdPP database
"""


import os
import errno
from sys import getsizeof as sz
import sqlite3
from datetime import datetime
from typing import Optional, List, Union, Tuple, Any, Iterator, Dict

import numpy as np
from numpy import typing as npt
from mzapy.isotopes import OrderedMolecularFormula
from rdkit import Chem
from rdkit import RDLogger

from idpp import __version__ as IDPP_VER
from idpp.msms.spectra import spec_combine


# TODO: Compose some complex type aliases to make some type annotations more understandable


# minimum Python version for proper DB support
_MIN_PYTHON_VER: str = "3.12.1"


# placeholder to use when no CE is provided for a MS/MS spectrum
_NO_CE_PLACEHOLDER: str = "_"


# TODO: It would be better to define a function that uses some patterns/rules
#       to do adduct normalization rather than using a statically defined set
#       of explicit replacements as is done here.

# statically define some adduct type replacements to apply in database insert
_REPLACE_ADDUCTS: Dict[str, str] = {
    # "value_to_replace": "replacement"
    "[M+CHO2]-": "[M+HCOO]-",
    "M-H": "[M-H]-",
    "M+H": "[M+H]+",
    "M+Na": "[M+Na]+",
    "[M]+*": "[M]+",
    "[M+Li]+*": "[M+Li]+",
    "M+K": "[M+K]+",
    "[M-H]1-": "[M-H]-",
    "[M+H]+[-H2O]": "[M+H-H2O]+",
    "[M-H2O+H]+": "[M+H-H2O]+",
    "[M-H2O-H]-": "[M-H-H2O]-",
    "[M-2H]-2": "[M-2H]2-",
    "[M+4H]+4": "[M+4H]4+",
    "[M+5H]+5": "[M+5H]5+",
    "[M+6H]+6": "[M+6H]6+",
    "[M+7H]+7": "[M+7H]7+",
    "[M+8H]+8": "[M+8H]8+",
    "[M+9H]+9": "[M+9H]9+",
    "[M+10H]+10": "[M+10H]10+",
    "[M+11H]+11": "[M+11H]11+",
    "[M+12H]+12": "[M+12H]12+",
    "[M+13H]+13": "[M+13H]13+",
    "[M+14H]+14": "[M+14H]14+",
    "[M+15H]+15": "[M+15H]15+",
    "[M-H]": "[M-H]-",
    "[M+2H]++": "[M+2H]2+",
    "[M+H+Na]+2": "[M+H+Na]2+",
    "[M+H+K]+2": "[M+H+K]2+",
    "M+": "[M]+",
    "[M]++": "[M]2+",
    "[M+]": "[M]+",
    "M+H-H2O": "[M+H-H2O]+",
    "M+2Na": "[M+2Na]2+",
    "M+Cl": "[M+Cl]-"
}


def _get_tstamp() -> str:
    """ returns a standardized timestamp (format: YY/MM/DD-hh:mm) """
    return datetime.now().strftime("%y/%m/%d-%H:%M")


def _db_ver_from_tstamp(tstamp: str
                        ) -> str :
    """
    convert a standard timestamp (i.e., from _get_tstamp()) into the
    corresponding database version string (format: YYMMDD.HH.MM)
    """
    return tstamp.replace("/", "").replace("-", ".").replace(":", ".")


def _add_version_info_and_change_log_entry(cur: sqlite3.Cursor
                                           ) -> None :
    """ 
    insert version info and change log entry into a new database using the provided `sqlite3.Cursor` 
    """
    qry_vi = "INSERT INTO VersionInfo VALUES (?,?,?);"
    cur.execute(qry_vi, (_MIN_PYTHON_VER, IDPP_VER, _db_ver_from_tstamp(_get_tstamp())))
    qry_cl = "INSERT INTO ChangeLog VALUES (?,?,?);"
    cur.execute(qry_cl, (_get_tstamp(), "idpp.db.util.create_db", "create database"))


def create_db(f: str, 
              overwrite: bool = False
              ) -> None :
    """
    creates a sqlite database for storing molecular properties

    raises a RuntimeError if the database already exists unless overwrite is True

    Parameters
    ----------
    f : ``str``
        filename/path of the database
    overwrite : ``bool``, default=False
        if the database file already exists and this flag is True, then overwrite existing database 
        and do not raise the RuntimeError
    """
    # see if the file exists
    if os.path.exists(f):
        if overwrite:
            os.remove(f)
        else:
            msg = "create_db: database file ({}) already exists"
            raise RuntimeError(msg.format(f))
    # initial connection creates the DB
    con = sqlite3.connect(f)  
    cur = con.cursor()
    # execute SQL script to set up the database
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../_include/idpp_db.sqlite3"), "r") as sql_f:
        cur.executescript(sql_f.read())
    # add the version information
    # and insert an entry into the change log for database creation
    _add_version_info_and_change_log_entry(cur)
    # save and close the database
    con.commit()
    con.close()


class IdPPdb():
    """
    Object for interacting with IdPP database

    Attributes
    ----------
    db_path : ``str``
        path to IdPP database file
    read_only : ``bool``
        database is read-only
    enforce_idpp_ver : ``bool``
            raise a RuntimeError if version of this package does not match idpp_ver in the database
    cur : ``sqlite3.Cursor``
        cursor for making queries directly from the underlying database
    uncommitted_changes : ``bool``
        flag indicating whether there are uncommitted changes to the database
    version_info : ``dict(str:str)``
        version information from the database (key: component, value: version)
    change_log : ``list(dict(str:str))``
        database change log (list of dicts with 'tstamp', 'author' and 'notes' entries)
    last_qry : ``str``
        stores the last query that was run, useful reference for the SQL queries
        that are being run behind the scenes and helpful for debugging
    check_insert_hits : ``int``
        count how many hits (i.e. an existing entry ID was returned) there were
        in this session from all insert method calls that use self._check_insert
    check_insert_misses : ``int``
        count how many misses (i.e. a new entry was added) there were
        in this session from all insert method calls that use self._check_insert
    last_check_insert_was_hit : ``bool``
        flag indicating if the last check_insert was a hit
    ledger_size : ``dict(str:int)``
        dict mapping the number of elements currently being stored (value) in the ledger 
        for each table (key)
    ledger_mem : ``int``
        estimated current memory footprint of the ledger in MB, this is a decent approximation
        of the total memory footprint of an IdPPdb instance since the ledger is by far the 
        largest component (assuming the interface was not initialized as read-only)
    combine_ms2 : ``bool``
        indicates whether MS2 spectra are combined on entry
    """

    # TODO: It might be technically more accurate to refer to "abundance" rather 
    #       than "intensity" when discussing MS/MS spectra. The "intensity" is 
    #       probably better suited to discussing the raw spectral data, and in
    #       this context working with normalized spectra "abundance" is a more 
    #       relevant term to use. This is not a critical point and these terms 
    #       do end up getting used somewhat interchangeably, but this change in
    #       terminology would make the descriptions slightly less ambiguous. This
    #       change would also need to be propagated across the package to avoid
    #       confusion, so it may end up being more work than it is worth.

    # TODO: Implement some fetch_X_by_id methods to fetch specific rows 
    #       from tables using their IDs  

    # TODO: Implement some variants of the fetch_X_data methods that return the
    #       data in like a polars dataframe or something like that. Might be easier
    #       to integrate into other people's workflows that way.

    # TODO: Is it possible/worth it to have spectrum combination be optional? 
    #           > Spectrum combination is now optional at the level of this interface
    #       Would this be more suited as some sort of configuration option somehow 
    #       stored within the database? 
    #           > Possibly, so I will keep these comments in place.

    # TODO: In the docstrings of all insert_X methods, indicate whether the method
    #       is a check_insert or nocheck_insert operation

    def __init__(self, 
                 db_path: str,
                 read_only: bool = False,
                 enforce_idpp_ver: bool = True,
                 combine_ms2: bool = False
                 ) -> None :
        """
        Create an instance of IdPPdb inteface object

        Parameters
        ----------
        db_path : ``str``
            path to IdPP database file
        read_only : ``bool``, default=False
            only allow fetch methods and do not maintain a rowid ledger,
            is faster to use if not modifying the database
        enforce_idpp_ver : ``bool``, default=True
            raise a RuntimeError if version of this package does not match idpp_ver in the database
        combine_ms2 : ``bool``, default=True
            combine MS2 spectra on inserting into the database
        """
        # store db path and flags
        self.__db_path = db_path
        self.__read_only = read_only
        self.__enforce_idpp_ver = enforce_idpp_ver
        #  connect to database
        self.__con = (
            sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True) if self.read_only 
            else sqlite3.connect(f"file:{self.db_path}?mode=rw", uri=True) 
        )
        self.__cur = self.__con.cursor()
        # keep track of whether there are uncommitted changes to the database
        self.__uncommitted_changes = False 
        # fetch version info and changelog
        self.__version_info = {} 
        self.__change_log = []
        self._fetch_version_info_and_change_log()  # sets self.__version_info and self.__change_log
        # last query starts empty
        self.__last_qry = ""
        # fill the ledger with existing rowids from specified tables 
        self.__ledger = {}
        if not self.__read_only:
            self._fill_ledger()
        # initialize check_insert counters to 0 for this session
        self.__check_insert_hits = 0
        self.__check_insert_misses = 0
        self.__last_check_insert_was_hit = False
        # store the combine_ms2 flag
        self.__combine_ms2 = combine_ms2

    def __repr__(self):
        """ string representation of this instance """
        s = (
            "IdPP("
            f"db_path='{self.__db_path}', "
            f"read_only={self.__read_only}, "
            f"enforce_idpp_ver={self.__enforce_idpp_ver},",
            f"combine_ms2={self.__combine_ms2}"
            ")"
        )
        return s

    def __str__(self):
        """ same as __repr__ for now """
        return self.__repr__()

    @property
    def db_path(self):
        return self.__db_path
    
    @property
    def read_only(self):
        return self.__read_only
    
    @property
    def enforce_idpp_ver(self):
        return self.__enforce_idpp_ver

    @property
    def cur(self):
        return self.__cur
    
    @property
    def uncommitted_changes(self):
        return self.__uncommitted_changes
    
    @property
    def version_info(self):
        return self.__version_info

    @property
    def change_log(self):
        return self.__change_log

    @property
    def last_qry(self):
        return self.__last_qry

    @property
    def check_insert_hits(self):
        return self.__check_insert_hits
    
    @property
    def check_insert_misses(self):
        return self.__check_insert_misses
    
    @property
    def last_check_insert_was_hit(self):
        return self.__last_check_insert_was_hit
    
    @property
    def ledger_size(self):
        return {k: len(v) for k, v in self.__ledger.items()}
    
    @property
    def ledger_mem(self):
        return self._estimate_ledger_size_mb()
    
    @property
    def combine_ms2(self): 
        return self.__combine_ms2
        
    def commit(self
               ) -> None :
        """ commit changes to the databse (write to file) """
        self.__con.commit()
        # unset the uncommitted changes flag if it had been set
        self.__uncommitted_changes = False

    def close(self,
              ignore_uncommitted_changes: bool = False,
              ) -> None :
        """ 
        close connection to the database, release the ledger memory (if not read-only) 

        Parameters
        ----------
        ignore_uncommitted_changes : ``bool``, default=False
            A RuntimeError is raised if the uncommitted changes flag is set (i.e., trying 
            to close the database without committing the changes that have been made), this
            option ignores that check and no error is raised
        """
        # NOTE: This check occurs before going through the actual closing steps. This way 
        #       it is possible to go back and call commit() method then try to close again.
        # handle uncommitted changes
        if not ignore_uncommitted_changes and self.__uncommitted_changes:
            msg = "IdPPdb: close: database connection closed without committing changes"
            raise RuntimeError(msg)
        # close DB connection
        self.__con.close()
        # get rid of the ledger (it can take up a lot of memory) just in case
        # an application using this interface has other stuff it needs to do 
        # after this interface instance has been closed but it may not have gone
        # out of scope yet and gotten GCed
        if not self.__read_only:
            del self.__ledger

    def vacuum(self
               ) -> None:
        """
        run the VACUUM command on the database, to shrink the file size if a lot of 
        data is deleted from it
        """
        self.__con.execute("VACUUM")

    def update_version_info(self
                            ) -> None :
        """
        update version information stored in the database to be consistent with the current
        version of this idpp package (and its expected minimum Python version)
        """
        # update database
        self.__cur.execute("DELETE FROM VersionInfo;")
        qry = "INSERT INTO VersionInfo VALUES (?,?,?);"
        self.__cur.execute(qry, 
                           (_MIN_PYTHON_VER, IDPP_VER, self.__version_info["db_ver"]))
        self.__last_qry = qry
        # set the uncommitted changes flag
        self.__uncommitted_changes = True
        # update this interface's instance vars
        self.__version_info["python_ver"] = _MIN_PYTHON_VER
        self.__version_info["idpp_ver"] = IDPP_VER

    def insert_change_log_entry(self,
                                author: str, 
                                notes: str
                                ) -> None:
        """
        insert an entry into the change log with automatically generated timestamp,
        also updates the db_ver column of the VersionInfo table

        Parameters
        ----------
        author : ``str``
            who is responsible for the changes
        notes : ``str``
            description of the changes
        """
        tstamp = _get_tstamp()
        # update the change log
        self.__cur.execute("INSERT INTO ChangeLog VALUES (?,?,?);", (tstamp, author, notes))
        self.__change_log.append({"tstamp": tstamp, "author": author, "notes": notes})
        # update db_ver
        db_ver = _db_ver_from_tstamp(tstamp)
        qry = "UPDATE VersionInfo SET db_ver=?"
        self.__cur.execute(qry, (db_ver,))
        self.__last_qry = qry
        self.__version_info["db_ver"] = db_ver
        # set the uncommitted changes flag
        self.__uncommitted_changes = True
    
    def insert_class_definition(self,
                                class_name: str, 
                                notes: Optional[str] = None
                                ) -> int:
        """
        insert an entry into the ClassDefs table, defining a classification

        Parameters
        ----------
        class_name : ``str``
            name for the classification (e.g. "lipid")
        notes : ``str``, optional
            additional description of the classification definition

        Returns
        -------
        cls_id : ``int``
            classification ID of the newly added class definition
        """
        return self._check_insert("ClassDefs", (class_name,), extra_vals=(notes,))

    def insert_class_label(self,
                           cls_id: int, cmpd_id: int
                           ) -> None :
        """
        insert a class label (specified by cls_id) for a compound (by cmpd_id)

        Parameters
        ----------
        cls_id : ``int``
            classification ID
        cmpd_id : ``int``
            compount ID to assign classification to 
        """
        # TODO: This might be better as a check_insert method
        _ = self._nocheck_insert("ClassLabels", (cls_id, cmpd_id), add_rowid_none=False)

    def insert_src(self, 
                   src_name: str, src_ref: str, 
                   src_notes: Optional[str] = None
                   ) -> int :
        """ 
        insert an entry into the Sources table, return the rowid
        (checks if it already exists first and returns existing rowid if so)
        
        Parameters
        ----------
        src_name : ``str``
            name for the source
        src_ref : ``str``
            link/reference info for the source
        src_notes : ``str``, optional
            optional notes/metadata for the source, in JSON format

        Returns
        -------
        src_id : ``int``
            identifier for the source that was just added
        """
        return self._check_insert('Sources', (src_name,), extra_vals=(src_ref, src_notes))
    
    def fetch_src_data(self, 
                       n_rows: int
                       ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from Sources table
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch

        Yields
        ------
        src_id : ``int``
            source identifier
        src_name : ``str``
            source name
        src_ref : ``str``
            source reference info
        src_notes : ``str``
            source notes (in JSON format)
        """
        yield from self._fetch_row_generator(('Sources',), ('src_id', 'src_name', 'src_ref', 'src_notes'), n_rows)

    def insert_form(self, 
                    form: str
                    ) -> int :
        """ 
        insert an entry into the Formulas table, return the rowid
        (checks if it already exists first and returns existing rowid if so)

        Always tries to convert to mzapy.isotopes.OrderedMolecularFormula first, 
        if that fails return -1 placeholder
        
        Parameters
        ----------
        form : ``str``
            molecular formula string

        Returns
        -------
        form_id : ``int``
            identifier for the formula that was just added (or already present)
        """
        try:
            # this produces molecular formulas with consistent atom ordering
            ord_form = str(OrderedMolecularFormula(form))
        except ValueError:
            # If there is some sort of problem with OrderedMolecualarFormula
            # parsing the original string, it will raise a ValueError. Handle
            # such cases by just returning the placeholder (-1) form_id 
            return -1
        except KeyError:
            # If the formula was parsable but contains an element that is not
            # recognized by OrderedMolecularFormula, the call to str() will cause
            # a KeyError. Handle this by proceeding with the unordered molecular 
            # formula, since this is an easier thing to fix in the future by patching 
            # the _ELEMENT_MONOISO_MASS constant in mzapy.isotopes
            ord_form = form
        return self._check_insert('Formulas', (ord_form,))
    
    def fetch_form_data(self, 
                        n_rows: int
                        ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from Formulas table
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch

        Yields
        ------
        form_id : ``int``
            formula identifier
        form : ``str``
            molecular formula
        """
        yield from self._fetch_row_generator(('Formulas',), ('form_id', 'form'), n_rows)

    def insert_smi(self, 
                   smi: str
                   ) -> int :
        """
        insert an entry into the Smiles table, return the rowid
        (checks if it already exists first and returns existing rowid if so)

        Parameters
        ----------
        smi : ``str``
            SMILES structure

        Returns
        -------
        smi_id : ``int``
            identifier for the SMILES structure that was just added (or already present)
        """
        # turn off rdkit logging messages
        RDLogger.DisableLog('rdApp.*') 
        if (mol := Chem.MolFromSmiles(smi)):
            # We can trust that if we were able to create a Mol object from the SMILES
            # structure, then there should be no problem generating a SMILES structure
            # from that Mol object and no need to check whether it worked.
            return self._check_insert("Smiles", (Chem.MolToSmiles(mol, 
                                                                  isomericSmiles=True, 
                                                                  canonical=True),))
        # fallback to placeholder ID for SMILES structures that did not work
        return -1
    
    def fetch_smi_data(self, 
                       n_rows):
        """
        fetch data from Smiles table
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch

        Yields
        ------
        smi_id : ``int``
            SMILES identifier
        smi : ``str``
            SMILES structure
        """
        yield from self._fetch_row_generator(('Smiles',), ('smi_id', 'smi'), n_rows)

    def insert_inchi(self, 
                     inchi_key: str, 
                     inchi: Optional[str] = None
                     ) -> int :
        """
        insert an entry into the InChIs table, return the rowid
        (checks if it already exists first and returns existing rowid if so)

        Parameters
        ----------
        inchi_key : ``str``
            InChI key
        inchi : ``str``, optional
            InChI structure string
        
        Returns
        -------
        inchi_id : ``int``
            identifier for the compound that was just added (or already present)
        """
        # TODO: If there is an entry for a given InChI key already present in the database
        #       but it does not have an associated full InChI structure and subsequently there
        #       is an attempt to add an entry with the same InChI key which DOES have an
        #       associated full InChI structure, then the new InChI structure ends up getting
        #       lost due to _check_insert only checking on the inchi_key. This method should 
        #       run the _check_insert first, then check to see if self.last_check_insert_was_hit
        #       is True. If so, then see if the corresponding existing entry has a full InChI
        #       structure and if not, then go ahead and update that with the new one. This would
        #       be helped out a lot by also implementing a fetch_inchi_data_by_id method, which
        #       can go and check existing entries from the _check_insert hits.
        # turn off rdkit logging messages
        RDLogger.DisableLog('rdApp.*') 
        # I am fairly sure that any molecule will only have a single InChI structure
        # so there is no need to regenerate them to ensure maximal collisions. The only
        # thing I am checking here is that the InChI key matches the InChI structure (if
        # provided) by regenerating the key and using that to add the entry
        if inchi is not None:
            if (ikey := Chem.InchiToInchiKey(inchi)):
                return self._check_insert("InChIs", (ikey,), extra_vals=(inchi,))
            else:
                # if we can't compute an InChI key then there was probably something 
                # wrong with the InChI structure
                return -1
        return self._check_insert("InChIs", (inchi_key,))
            
    def fetch_inchi_data(self, 
                         n_rows: int
                         ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from InChIs table
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch

        Yields
        ------
        inchi_id : ``int``
            InChI identifier
        inchi_key : ``str``
            InChI key
        inchi : ``str``
            full InChI structure
        """
        yield from self._fetch_row_generator(('InChIs',), ('inchi_id', 'inchi_key', 'inchi'), n_rows)
        
    def insert_cmpd(self, 
                    name: str, 
                    form_id: int = -1, smi_id: int = -1, inchi_id: int = -1
                    ) -> int :
        """
        insert an entry into the Compounds table, return the rowid
        (checks if it already exists first and returns existing rowid if so)
        required form_id, smi_id, and inchi_id default to placeholder values (-1)

        Parameters
        ----------
        name : ``str``
            compound name
        form_id : ``int``, default=-1
            Formulas row identifier
        smi_id : ``int``, default=-1
            Smiles row identifier
        inchi_id : ``int``, default=-1
            InChIs row identifier

        Returns
        -------
        cmpd_id : ``int``
            identifier for the compound that was just added (or already present)
        """
        return self._check_insert('Compounds', 
                                  (name,), 
                                  extra_vals=(form_id, smi_id, inchi_id),
                                  ignore_check_vals=[("",),])
    
    def fetch_cmpd_data(self, 
                        n_rows: int
                        ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from Compounds table (joined to Formulas, Smiles, and InChIs table)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch

        Yields
        ------
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        form_id : ``int``
            Formulas identifier
        form : ``str``
            molecular formula
        smi_id : ``int``
            Smiles identifier
        smi : ``str``
            SMILES structure
        inchi_id : ``int``
            InChIs identifier
        inchi_key : ``str``
            InChI key
        inchi : ``str``
            full InChI structure
        """
        tables = ('Compounds',
                  'JOIN Formulas USING(form_id)',
                  'JOIN Smiles USING(smi_id)',
                  'JOIN InChIs USING(inchi_id)')
        values = ('cmpd_id', 'cmpd_name', 
                  'form_id', 'form',
                  'smi_id', 'smi',
                  'inchi_id', 'inchi_key', 'inchi')
        yield from self._fetch_row_generator(tables, values, n_rows)
    
    def insert_adduct(self, 
                      adduct: str, cmpd_id: int, mz: float, z: int
                      ) -> int :
        """
        insert an entry into the Compounds table, return the rowid
        (checks if it already exists first and returns existing rowid if so)

        Parameters
        ----------
        adduct : ``str``
            adduct type
        cmpd_id : ``int``
            Compounds row identifier
        z : ``int``
            adduct charge
        mz : ``float``
            adduct m/z 
        
        Returns
        -------
        adduct_id : ``int``
            identifier for the adduct that was just added (or already present)
        """
        # TODO: Validate/standardize the format of adducts on insert
        if (replaced := _REPLACE_ADDUCTS.get(adduct)) is not None:
            adduct = replaced
        return self._check_insert('Adducts', (adduct, cmpd_id), extra_vals=(z, mz))
    
    def fetch_adduct_data(self, 
                          n_rows: int,
                          ionization: str = 'both'):
        """
        fetch data from Adducts table (joined to Compounds and Formulas table)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch
        ionization : ``str``, default='both'
            only yield adducts with a specified ionization state: +/pos/POS/etc. for positive, 
            -/NEG/neg/etc. for negative, or both (the default) for both

        Yields
        ------
        adduct_id : ``int``
            adduct identifier
        adduct : ``str``
            adduct
        adduct_z : ``int``
            adduct charge
        adduct_mz : ``real``
            adduct m/z
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        """
        # set the appropriate WHERE clause based on the ionization param
        ionization = self._parse_ionization(ionization)
        where = {
            '+': 'WHERE adduct_z > 0',
            '-': 'WHERE adduct_z < 0',
        }[ionization] if ionization in '+-' else ''
        tables = ('Adducts',
                  'JOIN Compounds USING(cmpd_id)')
        values = ('adduct_id', 'adduct',
                  'adduct_z', 'adduct_mz',
                  'cmpd_id', 'cmpd_name')
        yield from self._fetch_row_generator(tables, values, n_rows, where=where)

    def fetch_adduct_data_extended(self, 
                                   n_rows: int,
                                   ionization: str = "both", 
                                   require_smi: bool = False, 
                                   restrict_adducts: Optional[List[str]] = None
                                   ) -> Iterator[List[Tuple[Any]]] :
        """
        (extended) fetch data from Adducts table (joined to Compounds and Formulas table, 
        and additionally include all of the compound data: formula, smiles, inchi)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch
        ionization : ``str``, default='both'
            only yield adducts with a specified ionization state: +/pos/POS/etc. for positive, 
            -/NEG/neg/etc. for negative, or both (the default) for both
        require_smi : ``bool``, default=False
            add a constraint to the query to only return rows that have SMILES structures 
            (for the compound -> cmpd_smi)    
        restrict_adducts : ``list(str)``, optional
            if included, restrict the results to only include the specified adducts
        
        Yields
        ------
        adduct_id : ``int``
            adduct identifier
        adduct : ``str``
            adduct
        adduct_z : ``int``
            adduct charge
        adduct_mz : ``float``
            adduct m/z
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        cmpd_form_id : ``int``
            (compound) Formulas identifier
        cmpd_form : ``str``
            (compound) molecular formula
        cmpd_smi_id : ``int``
            (compound) Smiles identifier
        cmpd_smi : ``str``
            (compound) SMILES structure
        cmpd_inchi_id : ``int``
            (compound) InChIs id
        cmpd_inchi_key : ``str``
            (compound) InChI key
        cmpd_inchi : ``str``
            (compound) InChI structure
       """
        # set the appropriate WHERE clause based on the ionization param
        ionization = self._parse_ionization(ionization)
        where = {
            '+': 'WHERE adduct_z > 0',
            '-': 'WHERE adduct_z < 0',
        }[ionization] if ionization in "+-" else ""
        # optionally add a restriction to requre SMILES structures
        if require_smi:
            where = where + " AND " if where != "" else "WHERE "
            where += "smi IS NOT NULL"
        # optionally restrict the adducts to specified options
        if restrict_adducts is not None:
            adducts_str = ",".join([f"'{_}'" for _ in restrict_adducts])
            where = where + " AND " if where != "" else "WHERE "
            where += f"adduct IN ({adducts_str})"
        tables = ('Adducts',
                  'LEFT JOIN Compounds ON Adducts.cmpd_id=Compounds.cmpd_id', 
                  'LEFT JOIN Formulas AS cFormulas ON cFormulas.form_id=Compounds.form_id',
                  'LEFT JOIN Smiles USING(smi_id)',
                  'LEFT JOIN InChIs USING(inchi_id)')
        values = ('adduct_id', 'adduct',
                  'adduct_z', 'adduct_mz',
                  'Adducts.cmpd_id', 'cmpd_name',
                  'Compounds.form_id', 'cFormulas.form',
                  'smi_id', 'smi', 
                  'inchi_id', 'inchi_key', 'inchi')
        yield from self._fetch_row_generator(tables, values, n_rows, where=where)

    def fetch_adduct_id_by_cmpd_id(self,
                                   cmpd_id: int,
                                   restrict_adducts: Optional[List[str]] = None,
                                   ) -> Iterator[int] :
        """
        Fetch all adduct_ids associated with an input cmpd_id

        Parameters
        ----------
        cmpd_id : ``int``
            query compound ID

        Yields
        ------
        adduct_id : ``int``
            corresponding adduct ID
        """
        where = f"WHERE cmpd_id={cmpd_id} "
        if restrict_adducts is not None:
            where += f"AND adduct IN ({",".join(f'"{_}"' for _ in restrict_adducts)})"
        tables = ('Compounds',
                  'JOIN Adducts USING(cmpd_id)')
        values = ('adduct_id',)
        for row in self._fetch_row_generator(tables, values, -1, where=where):
            # unpack the query rows and just return the single int adduct_ids
            yield row[0]

    def insert_ccs(self, 
                   ccs: float, 
                   adduct_id: int, 
                   src_id: int,
                   ) -> int :
        """
        insert an entry into the CCSs table, return the rowid
        (DOES NOT check if already exists before adding, always adds new entry)

        Parameters
        ----------
        ccs : ``float``
            CCS value
        adduct_id : ``int``
            Adducts row identifier
        src_id : ``int``
            Sources row identifier

        Returns
        -------
        ccs_id : ``int``
            identifier for the CCS value that was just added
        """
        return self._nocheck_insert('CCSs', (ccs, adduct_id, src_id))
    
    def fetch_ccs_data(self, 
                       n_rows : int,
                       ionization: str = 'both', select_sources: Optional[List[Union[str, int]]] = None
                       ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from CCSs table (joined to Adducts, Compounds, and Sources tables)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch
        ionization : ``str``, default='both'
            only yield adducts with a specified ionization state: +/pos/POS/etc. for positive, 
            -/NEG/neg/etc. for negative, or both (the default) for both
        select_sources : ``list(str or int)``, optional
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        ccs_id : ``int``
            CCS identifier
        ccs : ``float``
            CCS value
        adduct_id : ``int``
            adduct identifier
        adduct : ``str``
            adduct
        adduct_z : ``int``
            adduct charge
        adduct_mz : ``real``
            adduct m/z
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        src_id : ``int``
            source identifier
        src_name : ``str``
            source name
        """
        # set the appropriate WHERE clause based on the ionization param
        ionization = self._parse_ionization(ionization)
        where = {
            '+': 'WHERE adduct_z > 0',
            '-': 'WHERE adduct_z < 0',
        }[ionization] if ionization in '+-' else ''
        # deal with selection of source IDs/names if select_sources was provided
        add_to_where = self._where_clause_from_select_sources(select_sources) if select_sources is not None else ''
        if where != '':
            where = where + ' AND ' + add_to_where if add_to_where != '' else where
        else:
            where = 'WHERE ' + add_to_where if add_to_where != '' else ''
        tables = ('CCSs',
                  'JOIN Adducts USING(adduct_id)',
                  'JOIN Compounds USING(cmpd_id)',
                  'JOIN Sources USING(src_id)')
        values = ('ccs_id', 'ccs',
                  'adduct_id', 'adduct', 'adduct_z', 'adduct_mz',
                  'cmpd_id', 'cmpd_name',
                  'src_id', 'src_name')
        yield from self._fetch_row_generator(tables, values, n_rows, where=where)

    def fetch_ccs_by_adduct_id(self, 
                               adduct_id : int,
                               select_sources: Optional[List[Union[str, int]]] = None
                               ) -> Iterator[Tuple[Any]] :
        """
        fetch data from CCSs table corresponding to a specified adduct_id
        yields one row at a time as tuples

        Parameters
        ----------
        adduct_id : ``int``
            query adduct ID
        select_sources : ``list(str or int)``, optional
            optionally restrict the sources for yielded rows
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        ccs_id : ``int``
            CCS identifier
        ccs : ``float``
            CCS value
        """
        where = f"WHERE adduct_id={adduct_id} "
        where += ("AND " + self._where_clause_from_select_sources(select_sources) 
                  if select_sources is not None 
                  else "")
        tables = ('CCSs',
                  'JOIN Sources USING(src_id)')
        values = ('ccs_id', 'ccs')
        yield from self._fetch_row_generator(tables, values, -1, where=where)

    def insert_rt(self, 
                  rt: float, 
                  adduct_id: int, 
                  src_id: int
                  ) -> int :
        """
        insert an entry into the RTs table, return the rowid
        (DOES NOT check if already exists before adding, always adds new entry)

        Parameters
        ----------
        rt : ``float``
            RT value
        adduct_id : ``int``
            Adducts row identifier
        src_id : ``int``
            Sources row identifier

        Returns
        -------
        rt_id : ``int``
            identifier for the RT value that was just added
        """
        return self._nocheck_insert('RTs', (rt, adduct_id, src_id))
    
    def fetch_rt_data(self, 
                      n_rows: int,
                      ionization: str = 'both', 
                      select_sources: Optional[List[Union[str, int]]] = None
                      ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from RTs table (joined to Adducts, Compounds, and Sources tables)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch
        ionization : ``str``, default='both'
            only yield adducts with a specified ionization state: +/pos/POS/etc. for positive, 
            -/NEG/neg/etc. for negative, or both (the default) for both
        select_sources : ``list(str/int)``, optional
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        rt_id : ``int``
            RT identifier
        rt : ``float``
            RT value
        adduct_id : ``int``
            adduct identifier
        adduct : ``str``
            adduct
        adduct_z : ``int``
            adduct charge
        adduct_mz : ``real``
            adduct m/z
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        src_id : ``int``
            source identifier
        src_name : ``str``
            source name
        """
        # set the appropriate WHERE clause based on the ionization param
        ionization = self._parse_ionization(ionization)
        where = {
            '+': 'WHERE adduct_z > 0',
            '-': 'WHERE adduct_z < 0',
        }[ionization] if ionization in '+-' else ''
        # deal with selection of source IDs/names if select_sources was provided
        add_to_where = self._where_clause_from_select_sources(select_sources) if select_sources is not None else ''
        if where != '':
            where = where + ' AND ' + add_to_where if add_to_where != '' else where
        else:
            where = 'WHERE ' + add_to_where if add_to_where != '' else ''
        tables = ('RTs',
                  'LEFT JOIN Adducts USING(adduct_id)',
                  'LEFT JOIN Compounds USING(cmpd_id)',
                  'LEFT JOIN Sources USING(src_id)')
        values = ('rt_id', 'rt',
                  'adduct_id', 'adduct', 'adduct_z', 'adduct_mz',
                  'cmpd_id', 'cmpd_name',
                  'src_id', 'src_name')
        yield from self._fetch_row_generator(tables, values, n_rows, where=where)

    def fetch_rt_by_cmpd_id(self, 
                            cmpd_id : int,
                            select_sources: Optional[List[Union[str, int]]] = None
                            ) -> Iterator[Tuple[Any]] :
        """
        fetch data from RTs table corresponding to a specified cmpd_id
        yields one row (rt_id, rt) at a time as tuples

        Parameters
        ----------
        cmpd_id : ``int``
            query compound ID
        select_sources : ``list(str or int)``, optional
            optionally restrict the sources for yielded rows
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        rt_id : ``int``
            RT identifier
        rt : ``float``
            RT value
        """
        where = f"WHERE cmpd_id={cmpd_id} "
        where += ("AND " + self._where_clause_from_select_sources(select_sources) 
                  if select_sources is not None 
                  else "")
        tables = ("RTs",
                  "JOIN Sources USING(src_id)",
                  "JOIN Adducts USING(adduct_id)")
        values = ("rt_id", "rt")
        yield from self._fetch_row_generator(tables, values, -1, where=where)

    def fetch_rt_by_adduct_id(self, 
                              adduct_id : int,
                              select_sources: Optional[List[Union[str, int]]] = None
                              ) -> Iterator[Tuple[Any]] :
        """
        fetch data from RTs table corresponding to a specified adduct_id
        yields one row (rt_id, rt) at a time as tuples

        Parameters
        ----------
        adduct_id : ``int``
            query adduct ID
        select_sources : ``list(str or int)``, optional
            optionally restrict the sources for yielded rows
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        rt_id : ``int``
            RT identifier
        rt : ``float``
            RT value
        """
        where = f"WHERE adduct_id={adduct_id} "
        where += ("AND " + self._where_clause_from_select_sources(select_sources) 
                  if select_sources is not None 
                  else "")
        tables = ("RTs",
                  "JOIN Sources USING(src_id)")
        values = ("rt_id", "rt")
        yield from self._fetch_row_generator(tables, values, -1, where=where)

    def insert_ms2(self, 
                   ms2_mz: npt.NDArray, 
                   ms2_i: npt.NDArray, 
                   adduct_id: int, 
                   src_id: int, 
                   ms2_ce: Optional[int] = None,
                   max_n_fragments: int = 256
                   ) -> int :
        """
        insert an entry into the MS2Spectra table, return the rowid

        If there is already a spectrum for the same adduct_id and self.combine_ms2 is set, 
        then combine the input spectrum with that one and update, otherwise add a new entry

        Parameters
        ----------
        ms2_mz : ``numpy.ndarray``
        ms2_i : ``numpy.ndarray``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of floats)
        adduct_id : ``int``
            Adducts row identifier
        src_id : ``int``
            Sources row identifier
        ms2_ce : ``int``, optional
            specify collision energy (voltage as int) for the spectrum
        max_n_fragments : ``int``, default=256
            only retain the top N fragments in each spectrum

        Returns
        -------
        ms2_id : ``int``
            identifier for the MS/MS spectrum that was just added
        """
        # TODO: log the queries
        # only bother checking for existing spectrum if self.combine_ms2 flag is set
        # check if a spectrum already exists then merge if it does
        if self.__combine_ms2 and (existing_spectrum := self._fetch_ms2_spectrum_by_adduct_id(adduct_id)) is not None:
            # get some metadata
            qry = "SELECT ms2_id, ms2_n_spectra, ms2_ce FROM MS2Spectra WHERE adduct_id=?"
            ms2_id, ms2_n_spectra, existing_ce = self.__cur.execute(qry, (adduct_id,)).fetchone()
            # merge existing spectrum with new spectrum
            ms2_mz, ms2_i = spec_combine([existing_spectrum, (ms2_mz, ms2_i)], [ms2_n_spectra, 1])
            # update spectrum metadata
            qry = "UPDATE MS2Spectra SET ms2_n_spectra=:n, ms2_ce=:c WHERE ms2_id=:i"
            # fallthrough: update CE string is still None
            self.__cur.execute(qry, {"n": ms2_n_spectra + 1, 
                                    "c": f"{existing_ce},{ms2_ce if ms2_ce is not None else _NO_CE_PLACEHOLDER}", 
                                    "i": ms2_id})
            # drop all fragments with that ms2_id
            qry = "DELETE FROM MS2Fragments WHERE ms2_id=?"
            self.__cur.execute(qry, (ms2_id,))
        else:
            # add a new entry to MS2Spectra with spectrum metadata
            qdata = (adduct_id, 1, str(ms2_ce) if ms2_ce is not None else _NO_CE_PLACEHOLDER)
            ms2_id = self._nocheck_insert("MS2Spectra", qdata)
        # limit the number of fragments per spectrum
        if len(ms2_mz) > max_n_fragments:
            idx = np.argsort(ms2_i)[::-1]
            ms2_mz = ms2_mz[idx][:max_n_fragments]
            ms2_i = ms2_i[idx][:max_n_fragments]
        # convert the (new or merged) spectrum into integer representation and add
        # all of the fragments to MS2Fragments
        for imz, ii in zip(*self._convert_spectrum_to_int_format(ms2_mz, ms2_i)):
            _ = self._nocheck_insert("MS2Fragments", (ms2_id, imz, ii), add_rowid_none=False)
        # add an MS2Sources entry
        # this should be check_insert
        _ = self._check_insert("MS2Sources", (ms2_id, src_id))
        # return the spectrum ID
        return ms2_id

    def fetch_ms2_data(self, 
                       n_rows: int,
                       ionization: str = 'both', 
                       select_sources: Optional[List[Union[str, int]]] = None
                       ) -> Iterator[List[Tuple[Any]]] :
        """
        fetch data from MS2Spectra, MS2Fragments, MS2Sources tables 
        (joined to Adducts and Compounds tables)
        yields one batch of rows at a time as list of tuples

        Parameters
        ----------
        n_rows : ``int``
            number of rows to yield in each batch
        ionization : ``str``, default='both'
            only yield adducts with a specified ionization state: +/pos/POS/etc. for positive, 
            -/NEG/neg/etc. for negative, or both (the default) for both
        select_sources : ``list(str/int)``, optional
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Yields
        ------
        ms2_id : ``int``
            MS2 spectrum identifier
        ms2_mz : ``numpy.ndarray(float)``
            m/z component of MS2 spectrum as an array
        ms2_i : ``numpy.ndarray(float)``
            intensity component of MS2 spectrum as an array
        ms2_ce : ``str``
            collision energy of the combined spectra
        adduct_id : ``int``
            adduct identifier
        adduct : ``str``
            adduct
        adduct_z : ``int``
            adduct charge
        adduct_mz : ``real``
            adduct m/z
        cmpd_id : ``int``
            compound identifier
        cmpd_name : ``str``
            compound name
        src_ids : ``list(int)``
            source identifiers
        """
        # set the appropriate WHERE clause based on the ionization param
        ionization = self._parse_ionization(ionization)
        where = {
            '+': 'WHERE adduct_z > 0',
            '-': 'WHERE adduct_z < 0',
        }[ionization] if ionization in '+-' else ''
        # deal with selection of source IDs/names if select_sources was provided
        add_to_where = self._where_clause_from_select_sources(select_sources) if select_sources is not None else ''
        if where != '':
            where = where + ' AND ' + add_to_where if add_to_where != '' else where
        else:
            where = 'WHERE ' + add_to_where if add_to_where != '' else ''
        # ?
        tables = ('MS2Spectra',
                  'LEFT JOIN Adducts USING(adduct_id)',
                  'LEFT JOIN Compounds USING(cmpd_id)',
                  'LEFT JOIN MS2Sources USING(ms2_id)',
                  'LEFT JOIN MS2Fragments USING(ms2_id)', 
                  'GROUP BY ms2_id')
        values = ('ms2_id', 
                  'GROUP_CONCAT(frag_imz)', 'GROUP_CONCAT(frag_ii)',
                  'ms2_ce',
                  'adduct_id', 'adduct', 'adduct_z', 'adduct_mz',
                  'cmpd_id', 'cmpd_name',
                  'GROUP_CONCAT(src_id)')
        for rows in self._fetch_row_generator(tables, values, n_rows, where=where):
            yield [row[:1] + 
                   self._convert_spectrum_from_int_format(
                        *[[int(_) for _ in row[i].split(",")] 
                          for i in [1, 2]],
                   ) + 
                   row[3:-1] + 
                   # this is a little icky but it's easier than writing a correct query
                   (list(set([int(_) for _ in row[-1].split(",")])),) 
                   for row in rows]
        
    def fetch_ms2_ids_by_adduct_id(self, 
                                   adduct_id : int,
                                   select_sources: Optional[List[Union[str, int]]] = None
                                   ) -> List[int] :
        """
        fetch IDs from MS2Spectra table corresponding to a specified adduct_id
        returns a list of the matching IDs

        Parameters
        ----------
        adduct_id : ``int``
            query adduct ID
        select_sources : ``list(str or int)``, optional
            optionally restrict the sources for yielded rows
            specify a list of specific sources (by source name, `str`, or sorce identifier, `int`) 
            to include, if `None` include all sources

        Returns
        -------
        ms2_ids : ``list(int)``
            list of matching ms2 IDs for the specified adduct ID
        """
        where = f"WHERE adduct_id={adduct_id} "
        where += ("AND " + self._where_clause_from_select_sources(select_sources) 
                  if select_sources is not None 
                  else "")
        tables = ("MS2Spectra",
                  "JOIN Sources USING(src_id)")
        values = ("ms2_id")
        return [
            _[0] for _ in self._fetch_row_generator(tables, values, -1, where=where)
        ]

    def insert_ext_id(self, 
                      cmpd_id: int, src_id: int, ext_id: str
                      ) -> None :
        """
        insert an entry into the ExternalIDs table

        *does not return a rowID since this table just associates external
        identifiers with existing compounds/sources*

        Parameters
        ----------
        cmpd_id : ``int``
            compound identifier
        src_id : ``int``
            source identifier
        ext_id : ``str``
            external identifier
        """
        _ = self._check_insert("ExternalIDs", (cmpd_id, src_id, ext_id))

    def insert_dataset(self, 
                       description: str, query: str
                       ) -> int :
        """
        insert an entry into the Datasets table

        Parameters
        ----------
        description : ``str``
            description of the dataset
        query : ``str``
            query used for selecting the dataset

        Returns
        -------
        dataset_id : ``int``
            identifier for the dataset that was just added
        """
        return self._nocheck_insert("Datasets", (description, query))
    
    def insert_analysis_result(self, 
                               dataset_id: int, 
                               counts: npt.NDArray,
                               mz_tol: float,
                               rt_tol: Optional[float] = None, 
                               ccs_tol: Optional[float] = None,
                               ms2_tol: Optional[float] = None
                               ) -> None :
        """
        Insert the results from one trial of probability analysis (i.e. the counts
        of matches for each compound using a specified set of search tolerances)
        into the AnalysisResults table of the database

        Parameters
        ----------
        dataset_id : ``int`` 
            dataset identifier
        counts : ``numpy.ndarray(numpy.int32)``
            counts of matches using specified tolerances for each compound in dataset
        mz_tol : ``float``
            m/z tolerance (in ppm)
        rt_tol : ``float``, optional
            RT tolerance (min), None if not used
        ccs_tol : ``float``, optional
            CCS tolerance (percent), None if not used
        ms2_tol : ``float``, optional
            MS/MS similarity threshold, None if not used
        """
        _ = self._nocheck_insert("AnalysisResults", 
                                 (dataset_id, len(counts), mz_tol, rt_tol, ccs_tol, ms2_tol, counts.tobytes()),
                                 add_rowid_none=False)

    def _release_and_maj_ver_match(self,
                                   pkg_idpp_ver: str,
                                   db_idpp_ver: str
                                   ) -> bool :
        """
        returns a boolean indicating if the package and database idpp versions share the same 
        release and major version

        when checking database version against interface version, only worry about release
        and major versions, ignore minor version. So `0.5.2` should be fine with `0.5.10` 
        but `0.5.2` and `0.6.0` are incompatible. Also additional version information (like
        for feature development versions, e.g. `0.6.12.dylan_0`) is ignored

        Parameters
        ----------
        pkg_idpp_ver : ``str``
        db_idpp_ver : ``str``
            package and database idpp versions ("release.major_version.minor_version")

        Returns
        -------
        match : ``bool``
            package and database idpp releases and major versions match
        """
        # should I check that both version strings are the correct format first?
        pkg_rel, pkg_maj, *_ = pkg_idpp_ver.split(".")
        db_rel, db_maj, *_ = db_idpp_ver.split(".")
        # only check releases and major versions, must be exact matches
        return pkg_rel == db_rel and pkg_maj == db_maj

    def _fetch_version_info_and_change_log(self, 
                                           ) -> None :
        """
        fetches the version information and change log from the database
        storing them in self.__version_info and self.__change_log

        raise a RuntimeError if the IdPP version of this package doesw not match what is in the database
        """
        # fetch and store version info
        py_ver, idpp_ver, db_ver = self.__cur.execute("SELECT python_ver, idpp_ver, db_ver FROM VersionInfo").fetchall()[0]
        self.__version_info["python_ver"] = py_ver
        self.__version_info["idpp_ver"] = idpp_ver
        self.__version_info["db_ver"] = db_ver
        # TODO: add some enforcement of minimum python version
        # enforce IdPP version of this interface matches what is in the database
        if self.__enforce_idpp_ver and not self._release_and_maj_ver_match(IDPP_VER, idpp_ver):
            # con is open at this point, close it first
            self.close()
            msg = ("IdPPdb: _fetch_version_info_and_change_log: "
                   f"the version of the installed idpp package is {IDPP_VER} "
                   f"but the idpp_ver in the database is: {idpp_ver}")
            raise RuntimeError(msg)
        # fetch all entries from the changelog
        for tstamp, author, notes in self.__cur.execute("SELECT tstamp, author, notes FROM ChangeLog"):
            self.__change_log.append({"tstamp": tstamp, "author": author, "notes": notes})
    
    def _fill_ledger(self
                     ) -> None :
        """
        fill the ledger with existing rowids for tables where they are tracked

        the ledger (self.__ledger) is a dict mapping table name to dicts mapping check
        values (tuples) to rowids (ints):

        .. code-block:: python3

            self.__ledger = {
                'Sources': {
                    # (src_name,): rowid
                    ('HMDB',): 1,
                    ('CCSbase',): 2,
                    ...
                },
                ...
            }
        
        """
        # TODO: In the future, keeping this complete ledger in memory might end up being kind
        #                    kind of a pain as the database continues to grow in size. It would probably be
        #                    beneficial to make this caching a bit more configurable so that a user could
        #                    choose to only load into memory the components that are needed for a given task.
        #                    For instance, if I am adding a dataset that I know does not contain any SMILES 
        #                    structures or InChIs, then there is no point in caching that data in-memory which
        #                    would save a lot of space since those are quite big. This would also necessitate 
        #                    implementing some fallback logic in _check_insert to deal with cases where there
        #                    is an attempt by a user to insert a value into an uncached table, but it would
        #                    probably still be pretty straight forward.

        # TODO: If the compound name has "|" in it, that denotes the presence of synonyms. These should be
        #       split out and represented separately in the ledger for efficient searching and check_insert
        #       will not need to know anything special about checking entries for inserting into Compounds.
        #       The insert_compounds method will need to be amended to accomodate splitting and searching
        #       for these synonyms.
        tables = {
            'Sources': {'sel_qry': 'SELECT src_id, src_name FROM Sources'},
            'Formulas': {'sel_qry': 'SELECT form_id, form FROM Formulas'},
            'Smiles': {'sel_qry': 'SELECT smi_id, smi FROM Smiles'},
            'InChIs': {'sel_qry': 'SELECT inchi_id, inchi FROM InChIs'},
            'Compounds': {'sel_qry': 'SELECT cmpd_id, cmpd_name FROM Compounds'},
            'Adducts': {'sel_qry': 'SELECT adduct_id, adduct, cmpd_id FROM Adducts'},
            "ClassDefs": {"sel_qry": "SELECT cls_id, cls_name, cls_desc FROM ClassDefs"},
            "ExternalIDs": {"sel_qry": "SELECT ROWID, cmpd_id, src_id, ext_id FROM ExternalIDs"},
            "MS2Sources": {"sel_qry": "SELECT ROWID, ms2_id, src_id FROM MS2Sources"}
        }
        for table, tbldata in tables.items():
            self.__ledger[table] = {}
            for rowid, *check_vals in self.__cur.execute(tbldata['sel_qry']):
                self.__ledger[table][tuple(check_vals)] = rowid
            # TODO: log these queries

    def _estimate_ledger_size_mb(self
                                 ) -> int :
        """
        estimate the current size (in MB) of the ledger, this is only a rough estimate but should be 
        good enough to get an idea for the memory footprint of the IdPPdb object since the ledger is 
        the biggest component

        Returns
        -------
        approx_ledger_size : ``int``
            approximate size of the ledger in MB
        """
        # size of the top-level dict itself and the dicts it contains
        # k is just a str so no need to go deeper
        # v is a dict so this component of the size only tracks the container size of the dicts
        # not their actual contents
        top_lvl_sz = sz(self.__ledger) + sum([sz(k) + sz(v) for k, v in self.__ledger.items()])
        # figure out the sizes of the actual contained items
        contained_sz = 0
        for table, d in self.__ledger.items():
            # the keys in d are tuples of ints/strings and the values are ints
            # contained size has 3 components:
            #   - size of tuple
            #   - size of tuple elements
            #   - size of value (int)
            contained_sz += sum([sz(k) + sum([sz(c) for c in k]) + sz(v) for k, v in d.items()])
        # sum together, convert to MB and return
        return int((top_lvl_sz + contained_sz) / (1024 * 1024))

    def _check_insert(self, 
                      table: str, 
                      check_vals: Tuple[Any], 
                      extra_vals: Optional[Tuple[Any]] = None,
                      ignore_check_vals: List[Tuple[Any]] = [],
                      ) -> int :
        """
        check if a value(s) already exists and to the
        specified table if not. Returns the ID from that table
        (either the existing rowid or the newly added rowid)

        Parameters
        ----------
        table : ``str``
            specify the table
        check_vals : ``tuple(...)``
            value(s) to check or insert
        extra_vals : ``tuple(...)``, optional
            additional values to add to the table along with the check_vals
        ignore_check_vals : ``list(tuple(...))``, default=[]
            Pass a list of check_vals tuples that should be ignored. This has the effect of always adding
            new entries for specific check_vals rather than only adding a single value and returning the
            existing rowid for subsequent adds. One use case for this behavior is compounds entries that 
            have empty strings for names (but that may have other associated identifiers or external 
            identifiers), which should be kept separate instead of all mapping to a single compound.

        Returns
        -------
        rowid : ``int``
            ID of the existing/newly added element
        """
        # construct the select and insert queries for the different tables
        # TODO: Describe the requirements for tbldata below, like what do all of the fields mean?
        tbldata = {
            'Sources': {'rowid': 'src_id', 'checkvals': ['src_name'], 'nvals': 4},
            'Formulas': {'rowid': 'form_id', 'checkvals': ['form'], 'nvals': 2},
            'Smiles': {'rowid': 'smi_id', 'checkvals': ['smi'], 'nvals': 2},
            'InChIs': {'rowid': 'inchi_id', 'checkvals': ['inchi_key'], 'nvals': 3},
            'Compounds': {'rowid': 'cmpd_id', 'checkvals': ['cmpd_name'], 'nvals': 5},
            'Adducts': {'rowid': 'adduct_id', 'checkvals': ['adduct', 'cmpd_id'], 'nvals': 5},
            "ClassDefs": {"rowid": "cls_id", "checkvals": ["cls_name"], "nvals": 3},
            "ExternalIDs": {"rowid": None, "checkvals": ["cmpd_id", "src_id", "ext_id"], "nvals": 3},
            "MS2Sources" : {"rowid": None, "checkvals": ["ms2_id", "src_id"], "nvals": 2}
        }[table]
        # check the ledger first
        rowid = self.__ledger[table].get(check_vals, None)
        self.__check_insert_hits += 1
        self.__last_check_insert_was_hit = True
        if rowid is None or check_vals in ignore_check_vals:
            # not in the database, add a new entry
            qry_ins = 'INSERT INTO {tbl} VALUES ({nvals});'.format(tbl=table, nvals=','.join(tbldata['nvals'] * '?'))
            # store query BEFORE executing, for debugging in case of an error or unexpected results
            self.__last_qry = qry_ins
            extra_vals = (None,) * (tbldata['nvals'] - 1 - len(tbldata['checkvals'])) if extra_vals is None else extra_vals
            qdata = check_vals + extra_vals
            if tbldata["rowid"] is not None:
                qdata = (None,) + qdata
            self.__cur.execute(qry_ins, qdata)
            rowid = self.__cur.lastrowid
            # update the ledger
            self.__ledger[table][check_vals] = rowid
            # undo increment self.__check_insert_hits and increment self.__check_insert_misses instead
            # doing it this silly way so no need to add an else: after this block just to increment
            # the correct counter
            self.__check_insert_hits -= 1
            self.__check_insert_misses += 1
            self.__last_check_insert_was_hit = False
            # set the uncommitted changes flag
            self.__uncommitted_changes = True
        return rowid
    
    def _nocheck_insert(self, 
                        table: str, vals: Tuple[Any],
                        add_rowid_none: bool = True
                        ) -> int :
        """
        adds values to specified table without checking if they exist first
        Returns newly added rowid

        Parameters
        ----------
        table : ``str``
            specify the table
        vals : ``tuple(...)``
            value(s) to insert
        add_rowid_none : ``bool``, default=True
            add a None at the front of the query data as a placeholder for rowid, so that 
            when the values are inserted the row is automatically assigned an identifier
            set to False for tables that do not have an autoincrementing identifier as 
            the first column

        Returns
        -------
        rowid : ``int``
            ID of the newly added element
        """
        # build insert query based on table
        tbldata = {
            "ClassDefs": {"nvals": 3},
            "ClassLabels": {"nvals": 2},
            "ExternalIDs": {"nvals": 3},
            "AdductsToSmiles": {"nvals": 2},
            "CCSs": {"nvals": 4},
            "RTs": {"nvals": 4},
            "MS2Spectra": {"nvals": 4},
            "MS2Fragments": {"nvals": 3},
            "Datasets": {"nvals": 3},
            "AnalysisResults": {"nvals": 7},
        }[table]
        qry_ins = "INSERT INTO {tbl} VALUES ({nvals});".format(tbl=table, nvals=",".join(tbldata["nvals"] * "?"))
        # store query BEFORE executing, for debugging in case of an error or unexpected results
        self.__last_qry = qry_ins
        # TODO: add_rowid_none could probably be directly included in the tbldata definition 
        #       above instead of being passed in as an argument
        if add_rowid_none:
            vals = (None,) + vals
        self.__cur.execute(qry_ins, vals)
        # set the uncommitted changes flag
        self.__uncommitted_changes = True
        return self.__cur.lastrowid
    
    def _fetch_row_generator(self, 
                             tables: Tuple[str], values: Tuple[str], n_rows: int,
                             where: str = ""
                             ) -> Iterator[Union[List[Tuple[Any]], Tuple[Any]]] :
        """
        create a generator that yields specified row values from a specified table
        one row at a time, used for fetch methods

        Parameters
        ----------
        tables : ``tuple(str)``
            table name (and optionally some JOIN clauses) to fetch rows from
        values : ``tuple(str)``
            values to fetch from the specified table
        n_rows : ``int``
            number of rows to yield in each batch, if -1 then yeild all one by one instead of 
            in batches
        where : ``str``, default=""
            additional clause for filtering the rows that are yielded
            'WHERE <condition(s)>' in SQL

        Yields
        ------
        rows : ``list(tuple(...))``
            batch of specified row values
        """
        qry = 'SELECT {vals} FROM {tbls} {whr};'.format(vals=",".join(values), 
                                                        tbls=" ".join(tables), 
                                                        whr=where)
        # store query BEFORE executing, for debugging in case of an error or unexpected results
        self.__last_qry = qry
        self.__cur.execute(qry)
        if n_rows == 0 or n_rows < -1:
            msg = f"IdPPdb: _fetch_row_generator: n_rows must be a positive number or -1 (was: {n_rows})"
            raise ValueError(msg)
        if n_rows == -1:
            while (row := self.__cur.fetchone()) is not None:
                yield row
        else:
            while (rows := self.__cur.fetchmany(size=n_rows)) != []:
                yield rows

    def _parse_ionization(self, 
                          ionization: str
                          ) -> str :
        """
        multiple methods use an ionization parameter for controlling rows that are selected, this
        can have the value of 'both' or a number of equivalent values for either positive or negative
        so to reduce duplication, this method parses those equivalent values and returns 'both', '+',
        or '-'. Raises an error if ionization cannot be parsed
        """
        if ionization == 'both':
            return ionization
        if len(ionization) == 1:
            if ionization in '+-':
                return ionization
        if len(ionization) >= 3:
            first3 = ionization[:3].lower()
            if first3 in ['pos', 'neg']:
                return {'pos':'+', 'neg':'-'}[first3]
        # if it didnt match any of the above cases, raise an error
        # con is open at this point, close it first
        self.close()
        raise ValueError(f"IdPPdb: _parse_ionization: ionization '{ionization}' not recognized") 

    def _where_clause_from_select_sources(self, 
                                          select_sources):
        """
        some fetch methods use a select_sources parameter which can be a list of ``int`` and/or
        ``str`` for selecting only specific source IDs or names, respectively. These translate
        to WHERE clauses and since the logic for doing that is a bit much to duplicate in all
        of those methods, this method takes care of it, returning a string with the appropriate
        WHERE clause based on the select_sources parameter.
        """
        src_ids, src_names = [], []
        for src in select_sources:
            if type(src) is int:
                src_ids.append(src)
            elif type(src) is str:
                # add quotes to names
                src_names.append(f'"{src}"')
            else: 
                # con is open at this point, close it first
                self.close()
                msg = ("IdPPdb: _where_clause_from_select_sources: sources can be selected using "
                       f"source ID (int) or source name (str), {src} is invalid (type: {type(src)})")
                raise ValueError(msg)
        # create the selection conditions for src IDs/names, combine if needed
        id_cond = (
                "src_id IN (" + ",".join([str(_) for _ in src_ids]) + ")" 
                if len(src_ids) > 0 
                else ""
        )
        name_cond = (
            "src_name IN (" + ",".join(src_names) + ")" 
            if len(src_names) > 0 
            else ""
        )
        if id_cond != "":
            if name_cond != "":
                return "(" + id_cond + " OR " + name_cond + ")"
            else:
                return id_cond
        else:
            return name_cond

    def _convert_spectrum_to_int_format(self,
                                        ms2_mz: npt.NDArray, ms2_i: npt.NDArray
                                        ) -> Tuple[List[int], List[int]] :
        """
        Convert an MS/MS spectrum into integer format

        - m/z is multiplied by 10^5 then converted to integer
        - intensity is normalized (such that complete spectrum intensities sum to 1) and
            converted to ppm (integer)

        Parameters
        ----------
        ms2_mz : ``numpy.ndarray``
        ms2_i : ``numpy.ndarray``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of floats)

        Returns
        -------
        ms2_imz : ``list(int)``
        ms2_ii : ``list(int)``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of integers)
        """
        ms2_imz = np.round((1e5 * ms2_mz)).astype(np.int32)
        ms2_ii = np.round(1e6 * (ms2_i / sum(ms2_i))).astype(np.int32)
        # only keep values that are at least 1 ppm relative abundance
        idx = ms2_ii > 0
        # convert to lists and return
        return ms2_imz[idx].tolist(), ms2_ii[idx].tolist()
    
    def _convert_spectrum_from_int_format(self,
                                          ms2_imz: List[int], ms2_ii: List[int]
                                          ) -> Tuple[npt.NDArray, npt.NDArray] :
        """
        Convert an MS/MS spectrum from integer format back to normal float arrays

        Parameters
        ----------
        ms2_imz : ``list(int)``
        ms2_ii : ``list(int)``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of integers)
        
        Returns
        -------
        ms2_mz : ``numpy.ndarray``
        ms2_i : ``numpy.ndarray``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of floats)
        """
        # convert to arrays of floats
        msms_mz, msms_i = np.array([ms2_imz, ms2_ii], dtype=np.float64)
        # scale back to original ranges
        return msms_mz / 1e5, msms_i / 1e6

    def _fetch_ms2_spectrum_by_adduct_id(self,
                                         adduct_id: int
                                         ) -> Optional[Tuple[npt.NDArray, npt.NDArray]] :
        """
        fetch an MS2 spectrum by adduct ID, covert from int representation and return it

        returns None if no existing spectrum was found for the specified adduct_id

        Parameters
        ----------
        adduct_id : ``int``
        
        Returns
        -------
        ms2_mz : ``numpy.ndarray``
        ms2_i : ``numpy.ndarray``
            m/z and intensity components of MS/MS spectrum (as numpy arrays, of floats) if 
            an existing spectrum was found otherwise None
        """
        # TODO: log query
        qry = "SELECT frag_imz, frag_ii FROM MS2Spectra JOIN MS2Fragments USING(ms2_id) WHERE adduct_id=?"
        if (res := self.__cur.execute(qry, (adduct_id,)).fetchall()) != []:
            # convert query results into lists of m/z and intensity
            msms_imz, msms_ii = np.array(res).T
            # covert from integer representation and return
            return self._convert_spectrum_from_int_format(msms_imz, msms_ii)
        return None


# TODO: implement a function that checks the integrity of a database file
#                    making sure all of the correct tables are present and (maybe optionally?)
#                    the version info is consistent with whatever version is being used to 
#                    check the database (i.e. invoking this function)
