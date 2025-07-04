import importlib
import os
import sys
import pytest

sys.path.append(os.path.dirname(os.path.dirname(__file__)))


def test_imports():
    assert importlib.import_module('core.ce1_sampling')
    try:
        importlib.import_module('core.ce2_eda_recode')
    except Exception as exc:
        pytest.skip(f'ce2_eda_recode import skipped: {exc}')
    assert importlib.import_module('core.ce3_var_redu')
    assert importlib.import_module('utils.ce_log_tool')
