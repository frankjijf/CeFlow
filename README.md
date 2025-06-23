# CeFlow

CeFlow is a data preparation pipeline written in Python. It provides a Streamlit
interface for sampling (`CE1`) and exploratory data analysis with recoding
(`CE2`). The project includes modules that can be used programmatically as well
as an example application using the Titanic dataset.

## Features

- **Sampling (CE1)**: splits data into modeling, validation and test sets with
  optional bootstrap resampling.
- **EDA and Recoding (CE2)**: univariate analysis, missing value imputation and
  transformation helpers, producing Python recode scripts and profiling tables.
- **Streamlit Interface**: an interactive UI (`CE_IDE.py`) to run the workflow
  without writing code.

## Repository Structure

```
Archive/          Previous versions of the Streamlit interface.
CE_PY_CODE/       Core modules and the main Streamlit application.
  ├─ CE1_Sampling.py      Sampling utilities.
  ├─ CE2_EDA_Recode.py    EDA and recoding utilities.
  ├─ CE_Log_Tool.py       Simple logging helper.
  ├─ CE_IDE.py            Streamlit interface.
  └─ Titanic_test/        Example dataset and variable lists.
```

The directory `streamlit_output/` (created at runtime) holds results and logs
when running the Streamlit app.

## Requirements

Python 3.8 or later is recommended. Install dependencies with

```bash
pip install -r requirements.txt
```

## Running the Streamlit Application

Execute the following command from the repository root:

```bash
streamlit run CE_PY_CODE/CE_IDE.py
```

The interface allows you to upload a dataset (CSV or Excel) or use the provided
Titanic sample located in `CE_PY_CODE/Titanic_test`. Follow the sidebar
navigation to configure variables, run sampling and perform recoding. Outputs and
logs are written to `streamlit_output/` by default.
## Keeping the Repository Clean

The folder `CE_PY_CODE/streamlit_output/` holds temporary outputs such as logs
and recoded datasets. These files are generated when running the Streamlit
application and should not be committed to Git. The `.gitignore` file at the
repository root already excludes this directory. If you create additional
runtime files, append their patterns to `.gitignore` to keep the history clean.


## Example Usage Without Streamlit

The modules can also be imported directly in your own scripts:

```python
from CE_PY_CODE.CE1_Sampling import CE_Sampling
from CE_PY_CODE.CE2_EDA_Recode import CE_EDA_Recode
```

Each function accepts a configuration dictionary; see the source code for
details.

## License

This project is provided as-is without an explicit license. Contact the authors
if you wish to reuse the code or datasets.
