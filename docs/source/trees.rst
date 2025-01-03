``idpp.probability.trees``
=======================================
Module with data structures that facilitate identification probability analysis.


.. code-block:: python

    # type aliases
    type PropertyTree = Union[MzTree, RtTree, CcsTree, Ms2Tree]
    # map cmpd_id (int) to matching cmpd_ids (set(int))
    type QueryResult = Dict[int, Set[int]]  
    # map adduct_id (int) to similarity contributions from other adducts
    # similarity contributions: data arrays for COO matrix
    type Similarities = Tuple[
        npt.NDArray[np.int32], 
        npt.NDArray[np.int32], 
        npt.NDArray[np.float32]
    ]



``MzTree``
---------------------------------------

.. autoclass:: idpp.probability.trees.MzTree

Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.probability.trees.MzTree.__init__

.. autofunction:: idpp.probability.trees.MzTree.query_all

.. autofunction:: idpp.probability.trees.MzTree.query_all_gen

.. autofunction:: idpp.probability.trees.MzTree.save

.. autofunction:: idpp.probability.trees.MzTree.load_attrs


``CcsTree``
---------------------------------------

.. autoclass:: idpp.probability.trees.CcsTree

Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.probability.trees.CcsTree.__init__

.. autofunction:: idpp.probability.trees.CcsTree.query_radius_single

.. autofunction:: idpp.probability.trees.CcsTree.query_radius

.. autofunction:: idpp.probability.trees.CcsTree.query_all

.. autofunction:: idpp.probability.trees.CcsTree.save

.. autofunction:: idpp.probability.trees.CcsTree.load_attrs


``RtTree``
---------------------------------------

.. autoclass:: idpp.probability.trees.RtTree

Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.probability.trees.RtTree.__init__

.. autofunction:: idpp.probability.trees.RtTree.query_radius_single

.. autofunction:: idpp.probability.trees.RtTree.query_radius

.. autofunction:: idpp.probability.trees.RtTree.query_all

.. autofunction:: idpp.probability.trees.RtTree.save

.. autofunction:: idpp.probability.trees.RtTree.load_attrs



``Ms2Tree``
---------------------------------------

.. autoclass:: idpp.probability.trees.Ms2Tree

Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: idpp.probability.trees.Ms2Tree.__init__

.. autofunction:: idpp.probability.trees.Ms2Tree.precompute_similarities

.. autofunction:: idpp.probability.trees.Ms2Tree.query_all

.. autofunction:: idpp.probability.trees.Ms2Tree.save

.. autofunction:: idpp.probability.trees.Ms2Tree.load_attrs


Utility
---------------------------------------

.. autoclass:: idpp.probability.trees.DatasetQueries

.. autofunction:: idpp.probability.trees.load_tree

.. autofunction:: idpp.probability.trees.construct_mz_tree

.. autofunction:: idpp.probability.trees.construct_ccs_tree

.. autofunction:: idpp.probability.trees.construct_rt_tree

.. autofunction:: idpp.probability.trees.construct_ms2_tree

.. autofunction:: idpp.probability.trees.construct_ms2_tree_for_adduct_ids

.. autofunction:: idpp.probability.trees.construct_property_trees
