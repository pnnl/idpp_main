# IdPP - Identification Probability and Precision
Code repository for work on characterizing the relationship between identification probability and precision in reference-free metabolomics. 

- Dylan Ross (dylan.ross@pnnl.gov)
- [doi?](TODO)
- [license](LICENSE.txt)
- [disclaimer](DISCLAIMER.txt)


## Notebooks and scripts
Jupyter notebooks and Python scripts with code for assembling and cleaning the IdPP database, in addition to performing identification probability analysis and some assorted benchmarks/tests.

> __NOTE:__ 
> Many of these scripts and notebooks rely on additional data that is not present in this repository, and thus
> it unlikely that they will be capable of being sucessfully run after simply downloading from here. 
> The primary purpose of these scripts and notebooks is to serve as documentation of the various analyses 
> performed for this work (though some may be functional or can be made so with varying levels of effort). 

### Composite MS/MS spectrum tests 
- [match tests with composite spectra](notebooks_and_scripts/composite_spectrum_tests/match_tests.ipynb)
- [plotting for match test results](notebooks_and_scripts/composite_spectrum_tests/plot_match_test_results.ipynb)
- [analysis of individual and combined spectra by PCA](notebooks_and_scripts/composite_spectrum_tests/pca_comparison.ipynb)

### MS/MS benchmarking
- [setup databases for benchmarking](notebooks_and_scripts/msms_float_int_benchmarks/create_test_databases.py)
- [benchmark script](notebooks_and_scripts/msms_float_int_benchmarks/benchmark.py)

### IdPP database setup and expansion
- [perform a fresh build of the IdPP database](notebooks_and_scripts/idpp_database/rebuild_idpp_db.py)
- [clean the database](notebooks_and_scripts/idpp_database/clean_db.py)
- [expand database with predicted values](notebooks_and_scripts/idpp_database/expand_db.py)

### Identification probability analysis
- [figure depicting different search tolerances and instrument capabilities](notebooks_and_scripts/probability_analysis/tolerance_figure.ipynb)
- [summarize database coverage and contents](notebooks_and_scripts/probability_analysis/database_summaries.ipynb)
- [identification probability analysis on original database](notebooks_and_scripts/probability_analysis/idprob_reference_only.ipynb)
- [identification probability analysis on expanded database](notebooks_and_scripts/probability_analysis/idprob_expanded.ipynb)
- [expanded database analysis with more caching of intermediate results](notebooks_and_scripts/probability_analysis/idprob_expanded_with_caching.ipynb)

## Documentation 
Documentation for selected components of the IdPP codebase (_i.e._ the `idpp` Python package)

> current idpp package version: 0.12.24

### IdPP database
- [builder](docs/idpp_db_builder.pdf)
- [interface](docs/idpp_db_interface.pdf)

### Utilities
- [MS2 spectrum similarity](docs/msms_utils.pdf)
- [identification probability tree data structures](docs/probability_trees.pdf)

