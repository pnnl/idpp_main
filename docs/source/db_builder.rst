``idpp.db.builder``
=======================================
Subpackage for populating the IdPP database with data from literature sources.


Overview
---------------------------------------
The build scripts in this subpackage each provide a function for adding a literature
data source into the database.


.. code-block:: python
        
    from idpp.db.builder.hmdb import add_hmdb_chunks_to_idppdb
    from idpp.db.builder.mona import add_mona_chunks_to_idppdb
    from idpp.db.builder.nist20 import add_nist20_msms_to_idppdb
    from idpp.db.builder.report import add_report_datasets_to_idppdb
    from idpp.db.builder.ccs_compendium import add_ccs_compendium_to_idppdb
    from idpp.db.builder.ccsbase import add_ccsbase_to_idppdb
    from idpp.db.builder.metlin_ccs import add_metlin_ccs_to_idppdb


There are also utility functions for taking some of the larger datasets and chunking
them into more managable sizes.


.. note:: 

    The literature datasets are not housed in this repository
    and must be downloaded from their respective sources.


Build functions
---------------------------------------


``idpp.db.builder.hmdb``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.hmdb.hmdb_chunks_exist

.. autofunction:: idpp.db.builder.hmdb.chunk_hmdb_xml

.. autofunction:: idpp.db.builder.hmdb.add_hmdb_chunks_to_idppdb

``idpp.db.builder.mona``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.mona.mona_chunks_exist

.. autofunction:: idpp.db.builder.mona.chunk_mona_json

.. autofunction:: idpp.db.builder.mona.add_mona_chunks_to_idppdb

``idpp.db.builder.nist20``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.nist20.add_nist20_msms_to_idppdb

``idpp.db.builder.report``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.report.add_report_datasets_to_idppdb

``idpp.db.builder.ccs_compendium``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.ccs_compendium.add_ccs_compendium_to_idppdb

``idpp.db.builder.ccsbase``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.ccsbase.add_ccsbase_to_idppdb

``idpp.db.builder.metlin_ccs``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.db.builder.metlin_ccs.add_metlin_ccs_to_idppdb
