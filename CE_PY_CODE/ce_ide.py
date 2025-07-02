

"""
Streamlit interface for CeFlow pipeline.
"""


import os
import sys
import pandas as pd
import numpy as np
import streamlit as st



# Add CE_PY_CODE to search path
BASE_DIR = os.path.dirname(__file__)
sys.path.append(os.path.join(BASE_DIR, "CE_PY_CODE"))



from ce1_sampling import CE_Sampling
from ce2_eda_recode import CE_EDA_Recode
from ce3_var_redu import CE_Var_Redu
from ce_log_tool import printto
from functools import partial



# Streamlit page configuration
st.set_page_config(page_title="CeFlow", layout="wide")



# Initialize session_state
state = st.session_state
if "step" not in state:
    state.step = 1


# Sidebar navigation labels
step_labels = [
    "1. Load The Dataset",
    "2. Variable Lists",
    "3. CE1 Sampling",
    "4. CE2 EDA & Var Recode",
    "5. CE3 Var Reduce",
]
# Render radio bound to session_state
sel = st.sidebar.radio("Step", step_labels, key="sidebar_sel")
# Update step based on sidebar selection
state.step = step_labels.index(sel) + 1

# Callback: button click updates sidebar_sel


def go_to(step):
    """Update sidebar selection and workflow step."""
    state.sidebar_sel = step_labels[step - 1]


st.title("CeFlow Streamlit IDE")

# -------- Step 1: Upload and Preview Data --------
if state.step == 1:
    st.header("Step 1: Upload and Preview Data", divider=True)
    uploaded = st.file_uploader("Upload data (csv/xlsx)", type=["csv", "xlsx", "xls"])
    path_input = st.text_input("Or input file path")
    df = None
    if uploaded:
        df = (
            pd.read_csv(uploaded)
            if uploaded.name.lower().endswith(".csv")
            else pd.read_excel(uploaded)
        )
    elif path_input:
        df = (
            pd.read_csv(path_input)
            if path_input.lower().endswith(".csv")
            else pd.read_excel(path_input)
        )
    if df is not None:
        state.df = df
        st.subheader("Data preview (first 100 rows)")
        st.dataframe(df.head(100))
        st.button("Next: Variable Setup", on_click=go_to, args=(2,))

# -------- Step 2: Variable Configuration --------
elif state.step == 2:
    st.header("Step 2: Variable Lists", divider=True)

    # ---------- 0) Data check ----------
    df = state.get("df")
    if df is None:
        st.warning("Please complete Step 1 first")
        st.stop()

    # ---------- 1) Data preview ----------
    st.subheader("Data preview (first 10 rows)", divider=True)
    st.dataframe(df.head(10))

    # ---------- 2) ID & dependent variable ----------
    st.subheader("ID & Dependent Variable", divider=True)
    st.selectbox("ID [Unique ID]", df.columns, index=0, key="id_col")
    st.selectbox(
        "DEP_VAR [Your dependent variable]", df.columns, index=1, key="dep_var"
    )
    st.pills(
        "BINARY_DV [Your dependent variable type]",
        ["Y", "N"],
        default="Y",
        key="binary_dv",
    )

    id_col, dep_var = state.id_col, state.dep_var
    all_vars = [c for c in df.columns if c not in (id_col, dep_var)]

    # ---------- 3) Initial automatic categorization ----------
    if "assigned" not in state:
        num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        auto_bin = [c for c in num_cols if df[c].nunique() == 2]
        auto_cont = [c for c in num_cols if c not in auto_bin]
        obj_cols = df.select_dtypes(include=["object", "category"]).columns.tolist()

        state.assigned = {
            "binary": auto_bin,
            "ordinal": [],  # auto-detection can be added if needed
            "nominal": obj_cols,
            "continuous": auto_cont,
        }

    # ---------- 4) Write assigned values into each widget's session_state key ----------
    cats = ["binary", "ordinal", "nominal", "continuous"]
    for cat in cats:
        widget_key = f"sel_{cat}"
        if widget_key not in state:
            state[widget_key] = state.assigned.get(cat, [])

    # ---------- 5) Define callback to keep variables globally unique ----------
    def handle_change(changed_cat: str):
        """When a category changes, remove new variables from the others"""
        new_vals = set(state[f"sel_{changed_cat}"])
        for other in cats:
            if other == changed_cat:
                continue
            other_key = f"sel_{other}"
            state[other_key] = [v for v in state[other_key] if v not in new_vals]

        # Sync back to assigned
        state.assigned = {c: state[f"sel_{c}"] for c in cats}

    # ---------- 7) Four-column multiselect ----------
    labels = {
        "binary": "BIN_VARS [Binary variables]",
        "ordinal": "ORD_VARS [Ordinal variables]",
        "nominal": "NOM_VARS [Nominal variables]",
        "continuous": "CONT_VARS [Continuous variables]",
    }
    cols = st.columns(4)

    for cat, col in zip(cats, cols):
        with col:
            st.multiselect(
                label=labels[cat],
                options=all_vars,  # full list
                key=f"sel_{cat}",
                on_change=partial(handle_change, cat),  # deduplicate in callback
            )
    state.id_col_val = state.id_col
    state.dep_var_val = state.dep_var
    state.binary_dv_val = state.binary_dv

    # ---------- 8) Next step ----------
    st.button("Next: CE1 Sampling", on_click=go_to, args=(3,))
# -------- Step 3: CE1 Sampling --------
elif state.step == 3:
    st.header("Step 3: CE1 Sampling", divider=True)
    df = state.get("df")
    if df is None:
        st.warning("Please complete the previous steps")
    else:
        assigned = state.get("assigned", {})
        cont_vars = assigned.get("continuous", [])
        nom_vars = assigned.get("nominal", [])
        bin_vars = assigned.get("binary", [])
        ord_vars = assigned.get("ordinal", [])
        id_col = state.get("id_col_val")
        dep_var = state.get("dep_var_val")
        binary_dv = state.get("binary_dv_val")
        left, right = st.columns([1, 2], border=True)
        with left:
            with st.form("ce1_form", border=False):
                st.subheader("General & Sampling Variables", divider=True)

                label_out_dir = """**PATH_OUTPUT**  
                :blue-background[Your output destination folder]"""
                out_dir = st.text_input(
                    label_out_dir, os.path.join(BASE_DIR, "streamlit_output")
                )

                label_split_portion = """**SPLIT_PORTION**  
                :blue-background[Your split portion for modeling portion]"""
                split_portion = st.number_input(label_split_portion, value=0.7)

                label_exclusion_if = """**EXCLUTION_IF**  
                :blue-background[Data exclusion condition; set an actual condition if needed, e.g., inds[dep_var].isna()]"""
                exclusion_if = st.text_input(
                    label_exclusion_if, value="inds[dep_var].isna()"
                )

                with st.expander(
                    "Optional Macro Variables: These all have defaults that can be used",
                    expanded=False,
                ):

                    label_bootstrap = """**BOOTSTRAP**  
                    :blue-background[Request oversampling, bootstrap if responders size is small (Y/N). Binary DV only]"""
                    bootstrap = st.selectbox(label_bootstrap, ["Y", "N"], index=0)

                    label_oversampled_rr = """**OVERSAMPLED_RR**  
                    :blue-background[Your desired oversample response rate on the modeling part of data, 0.1-0.5]"""
                    oversampled_rr = st.number_input(label_oversampled_rr, value=0.1)

                    label_min_num_resp = """**MIN_NUM_RESP**  
                    :blue-background[Minimum number of responders if sampling]"""
                    min_num_resp = st.number_input(label_min_num_resp, value=2000)

                    label_seed = """**SEED**  
                    :blue-background[Random selection seed if sampling]"""
                    seed = st.number_input(label_seed, value=123456)

                submitted = st.form_submit_button("Run CE1 Sampling")
            if submitted:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log1 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
                lst1 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")
                config1 = {
                    "cont_vars": cont_vars,
                    "nom_vars": nom_vars,
                    "bin_vars": bin_vars,
                    "ord_vars": ord_vars,
                    "path_output": out_dir,
                    "inds": df,
                    "id": id_col,
                    "dep_var": dep_var,
                    "binary_dv": binary_dv,
                    "weight": "",
                    "split_portion": split_portion,
                    "exclusion_if": "inds[dep_var].isna()",
                    "DS_present": "N",
                    "prefix": "R1_",
                    "keep_list": [],
                    "bootstrap": bootstrap,
                    "oversampled_rr": oversampled_rr,
                    "min_num_resp": min_num_resp,
                    "seed": int(seed),
                }
                state.config = config1
                with printto(log=log1, lst=lst1) as logger:
                    sampled_df, rate = CE_Sampling(config=config1, logger=logger)
                state.sampled_df = sampled_df
                state.ce1_rate = rate
                state.log1 = log1
                state.lst1 = lst1
        with right:
            if "sampled_df" in state:
                st.subheader("Sampling result", divider=True)
                st.write(f"Sample rate: {state.ce1_rate}")
                st.dataframe(state.sampled_df.head(100))
                st.download_button(
                    "Download sampled_df",
                    state.sampled_df.to_csv(index=False),
                    "CE2_sampled_df.csv",
                )
                with st.expander("View CE1 log"):
                    st.code(open(state.log1, encoding="utf-8").read())
                st.button("Next: CE2 Recoding", on_click=go_to, args=(4,))
            else:
                st.info("Fill parameters on the left and run CE1.")

# -------- Step 4: CE2 Recoding --------
elif state.step == 4:
    st.header("Step 4: CE2 Recoding", divider=True)
    sampled = state.get("sampled_df")
    if sampled is None:
        st.warning("Please complete CE1 sampling first")
    else:
        left, right = st.columns([1, 2], border=True)
        with left:
            with st.form("ce2_form", border=False):
                st.subheader("General & Sampling Variables", divider=True)

                label_out_dir = """**PATH_OUTPUT**  
                :blue-background[Your output destination folder]"""
                out_dir = st.text_input(
                    label_out_dir, os.path.join(BASE_DIR, "streamlit_output")
                )

                with st.expander(
                    "Options are ER for Equal Response method", expanded=False
                ):
                    label_profiling = """**PROFILING**  
                    :blue-background[Request profiling report on all variables (Y/N)]"""
                    profiling = st.pills(
                        label_profiling, ["Y", "N"], default="Y", key="profiling"
                    )

                    label_missrate = """**MISSRATE**  
                    :blue-background[Maximum missing rate allowed]"""
                    missrate = st.number_input(
                        label_missrate, value=0.75, key="missrate"
                    )

                    label_concrate = """**CONCRATE**  
                    :blue-background[Maximum amount of file that can be in a single value of an independent variable]"""
                    concrate = st.number_input(
                        label_concrate, value=0.9, key="concrate"
                    )

                    label_valcnt = """**VALCNT**  
                    :blue-background[Maximum number of unique values allowed. Set to 0 to allow any number]"""
                    valcnt = st.number_input(label_valcnt, value=50, key="valcnt")

                    label_minbinnc = """**MINBINNC**  
                    :blue-background[Minimum count in a bin to be usable. Set to 0 to use minbinpct only]"""
                    minbinnc = st.number_input(
                        label_minbinnc, value=500, key="minbinnc"
                    )

                    label_minbinnp = """**MINBINNP**  
                    :blue-background[Minimum percent of file in a bin to be usable. Set to 0 to use minbincnt only]"""
                    minbinnp = st.number_input(
                        label_minbinnp, value=0.05, key="minbinnp"
                    )

                    label_talpha = """**TALPHA**  
                    :blue-background[T-Test significance level for collapse of bins]"""
                    talpha = st.number_input(label_talpha, value=0.05, key="talpha")

                    label_bonfer = """**BONFER**  
                    :blue-background[Do Bonferroni adjustment for talpha? (Y/N)]"""
                    bonfer = st.pills(
                        label_bonfer, ["Y", "N"], default="N", key="bonfer"
                    )

                    label_nom_method = """**NOM_METHOD**  
                    :blue-background[Recoding method for nominal variables: Binary, Index or Mean]"""
                    nom_method = st.text_input(
                        label_nom_method, "INDEX", key="nom_method"
                    )

                    label_pvalue = """**PVALUE**  
                    :blue-background[P-value threshold to include variable]"""
                    pvalue = st.number_input(label_pvalue, value=0.05, key="pvalue")

                    label_min_size = """**MIN_SIZE**  
                    :blue-background[Minimum missing group size for Equal Response imputation]"""
                    min_size = st.number_input(
                        label_min_size, value=500, key="min_size"
                    )

                    label_num_category = """**NUM_CATEGORY**  
                    :blue-background[Maximum number of categories for profiling variables]"""
                    num_category = st.number_input(
                        label_num_category, value=10, key="num_category"
                    )

                    label_equal_dist = """**EQUAL_DIST**  
                    :blue-background[Use equal distance for dividing variables into groups for profiling (Y/N)]"""
                    equal_dist = st.pills(
                        label_equal_dist, ["Y", "N"], default="N", key="equal_dist"
                    )

                    label_p_lo = """**P_LO**  
                    :blue-background[Lower percentile for checking constant value]"""
                    p_lo = st.number_input(label_p_lo, value=1, key="p_lo")

                    label_p_hi = """**P_HI**  
                    :blue-background[Upper percentile for checking constant value]"""
                    p_hi = st.number_input(label_p_hi, value=99, key="p_hi")

                    label_impmethodC = """**IMPMETHODC**  
                    :blue-background[What method to use for missing imputation for continuous variables?]"""
                    impmethodC = st.text_input(
                        label_impmethodC, "median", key="impmethodC"
                    )

                    label_stdmethodC = """**STDMETHODC**  
                    :blue-background[Standardization options: any method allowed in proc stdize or NO to skip]"""
                    stdmethodC = st.text_input(
                        label_stdmethodC, "STD", key="stdmethodC"
                    )

                    label_cap_flrC = """**CAP_FLRC**  
                    :blue-background[What method to use for missing imputation for continuous variables?]"""
                    cap_flrC = st.pills(
                        label_cap_flrC, ["Y", "N"], default="Y", key="cap_flrC"
                    )

                    label_transformationC = """**TRANSFORMATIONC**  
                    :blue-background[Include transformed variables in evaluation (Y/N)]"""
                    transformationC = st.pills(
                        label_transformationC,
                        ["Y", "N"],
                        default="Y",
                        key="transformationC",
                    )

                    label_impmethodO = """**IMPMETHODO**  
                    :blue-background[What method to use for missing imputation for ordinal variables?]"""
                    impmethodO = st.text_input(
                        label_impmethodO, "mean", key="impmethodO"
                    )

                    label_stdmethodO = """**STDMETHODO**  
                    :blue-background[Standardization options: any method allowed in proc stdize or NO to skip]"""
                    stdmethodO = st.text_input(label_stdmethodO, "No", key="stdmethodO")

                    label_cap_flrO = """**CAP_FLRO**  
                    :blue-background[Do you want to do capping/flooring to handle outliers? (Y/N)]"""
                    cap_flrO = st.pills(
                        label_cap_flrO, ["Y", "N"], default="N", key="cap_flrO"
                    )

                    label_transformationO = """**TRANSFORMATIONO**  
                    :blue-background[Include transformed variables in evaluation? (Y/N)]"""
                    transformationO = st.pills(
                        label_transformationO,
                        ["Y", "N"],
                        default="N",
                        key="transformationO",
                    )

                submitted2 = st.form_submit_button("Run CE2")
            if submitted2:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log2 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
                lst2 = os.path.join(log_dir, "02_CE_EDA_Recode_LST_File.lst")
                config2 = state.get("config", {})
                config2.update(
                    {
                        "inds": sampled,
                        "profiling": profiling,
                        "missrate": missrate,
                        "concrate": concrate,
                        "valcnt": valcnt,
                        "minbinnc": minbinnc,
                        "minbinnp": minbinnp,
                        "talpha": talpha,
                        "bonfer": bonfer,
                        "nom_method": nom_method,
                        "pvalue": pvalue,
                        "min_size": min_size,
                        "num_category": num_category,
                        "equal_dist": equal_dist,
                        "p_lo": p_lo,
                        "p_hi": p_hi,
                        "impmethodC": impmethodC,
                        "stdmethodC": stdmethodC,
                        "cap_flrC": cap_flrC,
                        "transformationC": transformationC,
                        "impmethodO": impmethodO,
                        "stdmethodO": stdmethodO,
                        "cap_flrO": cap_flrO,
                        "transformationO": transformationO,
                        # "": ,
                    }
                )
                with printto(log=log2, lst=lst2) as logger:
                    rec, prof = CE_EDA_Recode(
                        indsn=sampled, config=config2, logger=logger
                    )
                state.config = config2
                state.recoded_df = rec
                state.profile_df = prof
                state.log2 = log2
                state.lst2 = lst2
        with right:
            if "recoded_df" in state:
                st.subheader("Recoding result", divider=True)
                st.dataframe(state.recoded_df.head(100))
                st.download_button(
                    "Download recoded data",
                    state.recoded_df.to_csv(index=False),
                    "CE2_Recoded.csv",
                )
                st.subheader("Profiling", divider=True)
                st.dataframe(state.profile_df)
                st.download_button(
                    "Download profiling report",
                    state.profile_df.to_csv(index=False),
                    "Profile.csv",
                )
                with st.expander("View CE2 log"):
                    st.code(open(state.log2, encoding="utf-8").read())
                st.button("Next: CE3 Var Reduce", on_click=go_to, args=(5,))
            else:
                st.info("Fill parameters on the left and run CE2.")

# -------- Step 5: CE3 Var Reduce --------
elif state.step == 5:
    st.header("Step 5: CE3 Var Reduce", divider=True)
    rec = state.get("recoded_df")
    if rec is None:
        st.warning("Please complete CE2 sampling first")
    else:
        left, right = st.columns([1, 2], border=True)
        with left:
            with st.form("ce3_form", border=False):
                st.subheader("Variable Reduction macro variables", divider=True)

                label_out_dir = """**PATH_OUTPUT**  
                :blue-background[Your output destination folder]"""
                out_dir = st.text_input(
                    label_out_dir, os.path.join(BASE_DIR, "streamlit_output")
                )

                label_fast_opt = """**FAST_OPT**  
                :blue-background[Fast option will turn off all tests except multivariate regression]"""
                fast_opt = st.pills(label_fast_opt, ["Y", "N"], default="N", key="fast_opt")

                with st.expander(
                    "Optional Macro Variables: These all have defaults that can be used",
                    expanded=False,
                ):
                    label_samplesize = """**SAMPLESIZE**  
                    :blue-background[Request profiling report on all variables (Y/N)]"""
                    samplesize = st.number_input(
                        label_samplesize, value=50000, key="samplesize"
                    )

                    label_redu_weight = """**REDU_WEIGHT**  
                    :blue-background[Use weights in variables reduction? (Y/N)]"""
                    redu_weight = st.pills(
                        label_redu_weight, ["Y", "N"], default="N", key="redu_weight"
                    )

                    label_sources = """**SOURCES**  
                    :blue-background[Minimum number of sources to be selected]"""
                    sources = st.number_input(label_sources, value=3, key="sources")

                    st.markdown("##### :rainbow[**Univariate regression option**]")
                    label_univ_reg_flg = """**UNIV_REG_FLG**  
                    :blue-background[Use univariate regression to choose variables? (Y/N)]"""
                    univ_reg_flg = st.pills(
                        label_univ_reg_flg, ["Y", "N"], default="Y", key="univ_reg_flg"
                    )
                    label_maxpuni = """**MAXPUNI**  
                    :blue-background[Maximum p value correlation for selecting via univariate regression]"""
                    maxpuni = st.number_input(label_maxpuni, value=0.0001, key="maxpuni")

                    st.markdown("##### :rainbow[**Correlation option**]")
                    label_correlation_flg = """**CORRELATION_FLG**  
                    :blue-background[Use correlation to choose variables? (Y/N)]"""
                    correlation_flg = st.pills(
                        label_correlation_flg,
                        ["Y", "N"],
                        default="Y",
                        key="correlation_flg",
                    )
                    label_corrcut = """**CORRCUT**  
                    :blue-background[Minimum correlation between independent variable and dependent variable]"""
                    corrcut = st.number_input(label_corrcut, value=0.01, key="corrcut")

                    st.markdown("##### :rainbow[**Principal components option**]")
                    label_principal_flg = """**PRINCIPAL_FLG**  
                    :blue-background[Use principal components to choose variables? (Y/N)]"""
                    principal_flg = st.pills(
                        label_principal_flg, ["Y", "N"], default="Y", key="principal_flg"
                    )
                    label_nprin = """**NPRIN**  
                    :blue-background[Number of principal components desired]"""
                    nprin = st.number_input(label_nprin, value=10, key="nprin")
                    label_minprin = """**MINPRIN**  
                    :blue-background[Minimum factor correlation for selecting via principal component]"""
                    minprin = st.number_input(label_minprin, value=0.5, key="minprin")

                    st.markdown("##### :rainbow[**Cluster option**]")
                    label_cluster_flg = """**CLUSTER_FLG**  
                    :blue-background[Use cluster analysis to choose variables? (Y/N)]"""
                    cluster_flg = st.pills(
                        label_cluster_flg, ["Y", "N"], default="Y", key="cluster_flg"
                    )
                    label_maxc = """**MAXC**  
                    :blue-background[Number of clusters desired]"""
                    maxc = st.number_input(label_maxc, value=20, key="maxc")
                    label_maxratio = """**MAXRATIO**  
                    :blue-background[Maximum R Squared ratio for selecting variables via clustering]"""
                    maxratio = st.number_input(label_maxratio, value=0.5, key="maxratio")

                    st.markdown("##### :rainbow[**Linear regression option**]")
                    label_regression_flg = """**REGRESSION_FLG**  
                    :blue-background[Use linear regression to choose variables? (Y/N)]"""
                    regression_flg = st.pills(
                        label_regression_flg, ["Y", "N"], default="Y", key="regression_flg"
                    )
                    label_alphareg = """**ALPHAREG**  
                    :blue-background[Alpha level for forward selection in linear regression]"""
                    alphareg = st.number_input(label_alphareg, value=0.05, key="alphareg")

                    st.markdown(
                        "##### :rainbow[**Logistic regression option - only applicable if binary dependent variable**]"
                    )
                    label_logistic_flg = """**LOGISTIC_FLG**  
                    :blue-background[Use logistic regression to choose variables? (Y/N)]"""
                    logistic_flg = st.pills(
                        label_logistic_flg, ["Y", "N"], default="Y", key="logistic_flg"
                    )
                    label_alphalog = """**ALPHALOG**  
                    :blue-background[Alpha level for forward selection in logistic regression]"""
                    alphalog = st.number_input(label_alphalog, value=0.05, key="alphalog")

                    st.markdown("##### :rainbow[**Information value option**]")
                    label_information_flg = """**INFORMATION_FLG**  
                    :blue-background[Use information value to choose variables? (Y/N)]"""
                    information_flg = st.pills(
                        label_information_flg,
                        ["Y", "N"],
                        default="Y",
                        key="information_flg",
                    )
                    label_decile = """**DECILE**  
                    :blue-background[Number of groups to use when calculate information values]"""
                    decile = st.number_input(label_decile, value=20, key="decile")
                    label_infvcut = """**INFVCUT**  
                    :blue-background[Minimum information value for selecting via information value]"""
                    infvcut = st.number_input(label_infvcut, value=0.01, key="infvcut")

                    st.markdown(
                        "##### :rainbow[**Maximum correlation between independent variables option**]"
                    )
                    label_ind_correlation_flg = """**IND_CORRELATION_FLG**  
                    :blue-background[Exclude variables with high correlation to others? (Y/N)]"""
                    ind_correlation_flg = st.pills(
                        label_ind_correlation_flg,
                        ["Y", "N"],
                        default="Y",
                        key="ind_correlation_flg",
                    )
                    label_maxcorr = """**MAXCORR**  
                    :blue-background[Maximum correlation allowed between independent variables]"""
                    maxcorr = st.number_input(label_maxcorr, value=0.7, key="maxcorr")

                    st.markdown(
                        "##### :rainbow[**Maximum correlation to dependent variable option**]"
                    )
                    label_ind_dv_corr_flg = """**IND_DV_CORR_FLG**  
                    :blue-background[Exclude variables with high correlation to dependent variable? (Y/N)]"""
                    ind_dv_corr_flg = st.pills(
                        label_ind_dv_corr_flg,
                        ["Y", "N"],
                        default="Y",
                        key="ind_dv_corr_flg",
                    )
                    label_max_dv_corr = """**MAX_DV_CORR**  
                    :blue-background[Maximum correlation allowed to dependent variable]"""
                    max_dv_corr = st.number_input(
                        label_max_dv_corr, value=0.7, key="max_dv_corr"
                    )

                submitted3 = st.form_submit_button("Run CE3")
            if submitted3:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log3 = os.path.join(log_dir, "03_CE_Var_Reduce_Log_File.log")
                lst3 = os.path.join(log_dir, "03_CE_Var_Reduce_LST_File.lst")
                config3 = state.get("config", {})
                config3.update(
                    {
                        "fast_opt": fast_opt,
                        "samplesize": samplesize,
                        "redu_weight": redu_weight,
                        "sources": sources,
                        "univ_reg_flg": univ_reg_flg,
                        "maxpuni": maxpuni,
                        "correlation_flg": correlation_flg,
                        "corrcut": corrcut,
                        "principal_flg": principal_flg,
                        "nprin": nprin,
                        "minprin": minprin,
                        "cluster_flg": cluster_flg,
                        "maxc": maxc,
                        "maxratio": maxratio,
                        "regression_flg": regression_flg,
                        "alphareg": alphareg,
                        "logistic_flg": logistic_flg,
                        "alphalog": alphalog,
                        "information_flg": information_flg,
                        "decile": decile,
                        "infvcut": infvcut,
                        "ind_correlation_flg": ind_correlation_flg,
                        "maxcorr": maxcorr,
                        "ind_dv_corr_flg": ind_dv_corr_flg,
                        "max_dv_corr": max_dv_corr,
                        "ind_correlation_flg": ind_correlation_flg,
                        "maxcorr": maxcorr,
                        "ind_dv_corr_flg": ind_dv_corr_flg,
                        "max_dv_corr": max_dv_corr,
                        "weight_is_freq": False,
                    }
                )
                with printto(log=log3, lst=lst3) as logger:
                    var_redu_df, varlist_redu = CE_Var_Redu(
                        df=rec, config=config3, logger=logger
                    )
                state.config = config3
                state.var_redu_df = var_redu_df
                state.varlist_redu = varlist_redu
                state.log3 = log3
                state.lst3 = lst3
        with right:
            if "var_redu_df" in state:
                st.subheader("var_redu_df", divider=True)
                st.dataframe(state.var_redu_df)
                st.subheader("varlist_redu", divider=True)
                st.dataframe(state.varlist_redu)
                with st.expander("CE3 Log"):
                    st.code(open(state.log3, encoding="utf-8").read())
            else:
                st.info("Fill parameters on the left and run CE3.")
