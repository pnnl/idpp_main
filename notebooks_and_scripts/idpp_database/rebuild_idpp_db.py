"""
    Do a fresh rebuild of the idpp database
"""


from idpp.db.util import create_db, IdPPdb
from idpp.db.builder.hmdb import add_hmdb_chunks_to_idppdb
from idpp.db.builder.mona import add_mona_chunks_to_idppdb
from idpp.db.builder.nist20 import add_nist20_msms_to_idppdb
from idpp.db.builder.report import add_report_datasets_to_idppdb
from idpp.db.builder.ccs_compendium import add_ccs_compendium_to_idppdb
from idpp.db.builder.ccsbase import add_ccsbase_to_idppdb
from idpp.db.builder.metlin_ccs import add_metlin_ccs_to_idppdb


def _main():
    dbf = "idpp.db"
    create_db(dbf, overwrite=True)
    db = IdPPdb(dbf)
    add_hmdb_chunks_to_idppdb(db, "hmdb_chunks/")
    add_mona_chunks_to_idppdb(db, "mona_chunks/")
    add_nist20_msms_to_idppdb(db, "msms_2020.db")
    add_report_datasets_to_idppdb(db, "RepoRT-master/")
    add_ccs_compendium_to_idppdb(db, "UnifiedCCSCompendium_FullDataSet_2024-04-19.csv")
    add_ccsbase_to_idppdb(db, "C3S.db")
    add_metlin_ccs_to_idppdb(db, "METLIN-CCS-03-15-2024.xlsx")
    db.commit()
    db.close()


if __name__ == "__main__":
    _main()
