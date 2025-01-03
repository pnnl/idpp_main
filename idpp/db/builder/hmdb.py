"""
    idpp/db/builder/hmdb.py

    Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for adding data downloaded from HMDB (https://hmdb.ca/downloads)
    - Metabolite and Protein Data (in XML Format)
        - Data Set: All Metabolites
        - Released On: 2021-11-17
"""


import errno
import os
import glob
from xml.etree import ElementTree

from idpp.db.util import IdPPdb


def hmdb_chunks_exist(chunk_dir: str
                      ) -> bool :
    """
    Check that chunked HMDB download files exist
    """
    # if the directory exists and there is at least 1
    # of ####.xml files in there, then call it good
    return os.path.isdir(chunk_dir) and len(glob.glob("[0-9][0-9][0-9][0-9].xml", root_dir=chunk_dir)) > 0


def chunk_hmdb_xml(hmdb_metabolites_xml: str,
                   chunk_dir: str,
                   ) -> None :
    """
    process a download of "hmdb_metabolites.xml" from HMDB into
    chunks of a managable size (~5 MB)

    Parameters
    ----------
    hmdb_metabolites : ``str``
        path to "hmdb_metabolites.xml" (download from HMDB)
    chunk_dir : ``str``
        directory to save the chunked files into
    """
    # ensure that the hmdb_metabolites.xml download file exists
    if not os.path.isfile(hmdb_metabolites_xml):
        # TODO: Offer to download the file from HMDB? It is pretty big (a few GB) but
        #       could be pretty convenient to have that option.
        raise FileNotFoundError(errno.ENOENT, 
                                os.strerror(errno.ENOENT), 
                                hmdb_metabolites_xml)
    # ensure that the chunk directory exists, if not make it
    if not os.path.isdir(chunk_dir):
        os.mkdir(chunk_dir)
    update_freq = 5
    # break into ~5MB chunks 
    chunk_size = int(1024 * 1024 * 5)
    with open(hmdb_metabolites_xml, 'r', encoding="utf8") as f:
        # find the end of the file and store byte number
        end = f.seek(0, os.SEEK_END)
        # jump back to the beginning and scan to the first line containing <metabolite> tag
        f.seek(0)
        line = ''
        while line.strip() != '<metabolite>':
            line = f.readline()
        pos = f.tell()
        # at the start of each loop line = '<metabolite>\n'
        chunk = 1
        while pos < end:
            buf = '<chunk>\n'
            s = 0
            while s < chunk_size and line.strip() != '</hmdb>':  
                buf += line
                line = f.readline()
                s += len(line)
            # either line is '</hmdb>' which is end of file
            # or some random stuff in the middle of a record
            # in which case need to continue until the closing
            # </metabolite> tag
            if line.strip() != '</hmdb>':
                # continue getting lines until </metabolite> tag
                while line.strip() != '</metabolite>':
                    buf += line
                    line = f.readline()
                buf += line
                line = f.readline()
            buf += '</chunk>\n'
            # write chunk to file
            with open(os.path.join(chunk_dir, f"{chunk:04d}.xml"), 'w', encoding="utf8") as out:
                out.write(buf)
            chunk += 1
            # update position
            pos = f.tell()
            if chunk % update_freq == 0:
                print('\r{:5d} chunks written {}  '.format(chunk, '|/--\\'[chunk // update_freq % 5]), end='')
    print('\r{:5d} chunks written {}  '.format(chunk, '|/--\\'[chunk // update_freq % 5]), end='')
    print()
            

def add_hmdb_chunks_to_idppdb(db: IdPPdb,
                              chunk_dir: str
                              ) -> None :
    """
    Add chunked HMDB metabolite data to IdPPdb

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface 
    chunk_dir : ``str``
        directory containing chunked HMDB download XML files
    """
    print("Adding HMDB compounds to IdPPdb ...")
    # make sure chunk_dir is valid and has chunk files in it
    if not hmdb_chunks_exist(chunk_dir):
        msg = f"add_hmdb_chunks_to_idppdb: HMDB chunks not found in {chunk_dir}"
        raise RuntimeError(msg)
    # add sources: HMDB and PubChem
    # source notes NULL for now, can update later
    hmdb_src_id = db.insert_src('HMDB', 'https://hmdb.ca')  
    pchm_src_id = db.insert_src('PubChem', 'https://pubchem.ncbi.nlm.nih.gov/')
    i = 0
    for chunk in sorted(glob.glob("[0-9][0-9][0-9][0-9].xml", root_dir=chunk_dir)):
        # parse the XML document
        root = ElementTree.parse(os.path.join(chunk_dir, chunk)).getroot()
        for metabolite in root:
            # extract relevant info from each metabolite
            name = metabolite.find('name').text
            hmdb_id = metabolite.find('accession').text
            smi = metabolite.find('smiles').text
            inchi = metabolite.find('inchi').text
            inchi_key = metabolite.find('inchikey').text
            e_pcid = metabolite.find('pubchem_compound_id')
            pubchem_cid = e_pcid.text if e_pcid is not None else None
            form = metabolite.find('chemical_formula').text
            # TODO: grab classification info from taxonomy section
            class_labels = [f"{level}:{res.text}" 
                            for level in ["kingdom", "superclass", "class", "subclass", "direct parent"] 
                            if (res := metabolite.find('taxonomy').find(level)) is not None]
            # TODO: grab synonyms from synonyms section
            # add data to database 
            form_id = db.insert_form(form)
            smi_id = db.insert_smi(smi) if smi is not None else -1
            inchi_id = db.insert_inchi(inchi_key, inchi=inchi) if inchi_key is not None else -1
            cmpd_id = db.insert_cmpd(name, form_id, smi_id=smi_id, inchi_id=inchi_id)
            # add in external compound identifiers (where available)
            db.insert_ext_id(cmpd_id, hmdb_src_id, hmdb_id)
            if pubchem_cid is not None:
                db.insert_ext_id(cmpd_id, pchm_src_id, pubchem_cid)
            # add in classification info (if found)
            for class_label in class_labels:
                cls_id = db.insert_class_definition(class_label)
                db.insert_class_label(cls_id, cmpd_id)
            # report progress
            i += 1
        print('\r\t{:8d} entries processed (last chunk: {:8s}) '.format(i, chunk), end='')
    print('\r\t{:8d} entries processed'.format(i), end='')
    print()
    # add a change log entry
    db.insert_change_log_entry("idpp.db.builder.hmdb.add_hmdb_chunks_to_idppdb", 
                               "add HMDB compounds")
    print("... done")
