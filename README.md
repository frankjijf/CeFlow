# SAS2PYTHON

This repository contains a Python translation of the SAS CE macros used for data sampling,
exploratory data analysis and recoding.  The example implementation works on the Titanic
passenger dataset located under `CE_PY_CODE/Titanic_test`.

## Requirements

Install the required packages with:

```bash
pip install pandas numpy scipy statsmodels scikit-learn xlsxwriter
```

Python 3.8 or newer is recommended.

## Basic usage

Run the main script from the repository root:

```bash
python CE_PY_CODE/CE0_fun_call_v5.py
```

The script loads the dataset, performs sampling via `CE_Sampling`, then calls
`CE_EDA_Recode` to create recoded variables and a profiling report.  Results and logs are
written to the `ce_output` and `log` subdirectories under the configured data folder.

### Configuration

Edit `CE_PY_CODE/CE0_fun_call_v5.py` to point to your own data directory and adjust variable
lists if necessary.

```python
lib_dir = r"D:\OneDrive\CE_PROJECT\Python_CE_Project_GH\SAS2PYTHON\CE_PY_CODE\Titanic_test"
out_dir = os.path.join(lib_dir, "ce_output")

# Variable lists
"cont_vars": contvar_list.contvar_list(),
"nom_vars": nomvar_list.nomvar_list(),
"bin_vars": binvar_list.binvar_list(),
"ord_vars": ordvar_list.ordvar_list(),
```

The modules in `CE_PY_CODE/Titanic_test` (`contvar_list.py`, `nomvar_list.py`,
`binvar_list.py`, `ordvar_list.py`) define the variables used in the pipeline.  Modify these
or create your own modules to match your dataset.

After configuring the paths and variable lists, run the command above and the generated
reports will appear in your output directory.
