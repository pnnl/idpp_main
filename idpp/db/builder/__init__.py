"""
    idpp/db/builder/__init__.py

    Dylan Ross (dylan.ross@pnnl.gov)

    subpackage for database building
"""

# TODO: Add a module with a function that selects out groups of MS2 spectra with the 
#       same adduct_id then combines them and fixes all of the corresponding database 
#       entries with updated counts/sources/etc. and drops all the duplicates. This 
#       way, the database can be efficiently filled with initial (individual) spectra
#       without combining, then groups of spectra can be combined in one shot afterwards. 
#       Combining spectra on insert makes sense for adding small numbers of spectra but 
#       it does not scale well with large collections. This would make a nice workaround
#       for that problem.
