import importlib
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(__file__)))


def test_imports():
    assert importlib.import_module('CE_PY_CODE.ce1_sampling')
    assert importlib.import_module('CE_PY_CODE.ce2_eda_recode')
    assert importlib.import_module('CE_PY_CODE.ce_log_tool')
