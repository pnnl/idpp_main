"""
    idpp/test/__include.py

    Dylan Ross (dylan.ross@pnnl.gov)

    special module defining a constant that other test modules can use
    to refer to built-in non-Python files used for testing purposes
"""


import os


# store path to test subpackage _include dir
TEST_INCLUDE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "_include"))
