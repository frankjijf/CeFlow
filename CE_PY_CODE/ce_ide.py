"""Streamlit interface for CeFlow pipeline."""

import os
import sys
import pandas as pd
import numpy as np
import streamlit as st

# 将 CE_PY_CODE 目录加入搜索路径
BASE_DIR = os.path.dirname(__file__)
sys.path.append(os.path.join(BASE_DIR, "CE_PY_CODE"))

from ce1_sampling import CE_Sampling
from ce2_eda_recode import CE_EDA_Recode
from ce_log_tool import printto
from functools import partial

# Streamlit 页面配置
st.set_page_config(page_title="CeFlow", layout="wide")

# 初始化 session_state
state = st.session_state
if "step" not in state:
    state.step = 1

# 侧边栏导航 标签
step_labels = ["1. Load the dataset", "2. Variable Lists", "3. CE1", "4. CE2"]
# 渲染 radio，key 与 session_state 绑定，无需 index
sel = st.sidebar.radio("Step", step_labels, key="sidebar_sel")
# 根据侧边栏选择更新 step
state.step = step_labels.index(sel) + 1

# 回调：按钮点击更新 sidebar_sel


def go_to(step):
    """Update sidebar selection and workflow step."""
    state.sidebar_sel = step_labels[step - 1]


st.title("CeFlow Streamlit IDE")

# -------- Step 1: 上传和预览数据 --------
if state.step == 1:
    st.header("Step 1: 上传和预览数据", divider=True)
    uploaded = st.file_uploader("上传数据 (csv/xlsx)", type=["csv", "xlsx", "xls"])
    path_input = st.text_input("或输入文件路径")
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
        st.subheader("数据预览 (仅前 100 行)")
        st.dataframe(df.head(100))
        st.button("下一步: 变量确认", on_click=go_to, args=(2,))

# -------- Step 2: 变量确认与配置 --------
elif state.step == 2:
    st.header("Step 2: Variable Lists", divider=True)

    # ---------- 0) 数据检查 ----------
    df = state.get("df")
    if df is None:
        st.warning("请先完成 Step 1 上传数据")
        st.stop()

    # ---------- 1) 数据预览 ----------
    st.subheader("数据预览 (前 10 行)", divider=True)
    st.dataframe(df.head(10))

    # ---------- 2) 标识 & 因变量 ----------
    st.subheader("标识 & 因变量配置", divider=True)
    st.selectbox("ID [Unique ID]", df.columns, index=0, key="id_col")
    st.selectbox("DEP_VAR [Your dependent variable]", df.columns, index=1, key="dep_var")
    st.pills("BINARY_DV [Your dependent variable type]", ["Y", "N"], default="Y", key="binary_dv")

    id_col, dep_var = state.id_col, state.dep_var
    all_vars = [c for c in df.columns if c not in (id_col, dep_var)]

    # ---------- 3) 初次自动分类 ----------
    if "assigned" not in state:
        num_cols  = df.select_dtypes(include=[np.number]).columns.tolist()
        auto_bin  = [c for c in num_cols if df[c].nunique() == 2]
        auto_cont = [c for c in num_cols if c not in auto_bin]
        obj_cols  = df.select_dtypes(include=["object", "category"]).columns.tolist()

        state.assigned = {
            "binary":      auto_bin,
            "ordinal":     [],          # 需要时可补自动识别
            "nominal":     obj_cols,
            "continuous":  auto_cont,
        }

    # ---------- 4) 把 assigned 值写进每个 widget 的 session_state 键 ----------
    cats = ["binary", "ordinal", "nominal", "continuous"]
    for cat in cats:
        widget_key = f"sel_{cat}"
        if widget_key not in state:
            state[widget_key] = state.assigned.get(cat, [])

    # ---------- 5) 定义回调，使变量全局唯一 ----------
    def handle_change(changed_cat: str):
        """当某类别 multiselect 变动时，把新增变量从其他类别里剔除"""
        new_vals = set(state[f"sel_{changed_cat}"])
        for other in cats:
            if other == changed_cat:
                continue
            other_key = f"sel_{other}"
            state[other_key] = [v for v in state[other_key] if v not in new_vals]

        # 同步回 assigned
        state.assigned = {c: state[f"sel_{c}"] for c in cats}

    # ---------- 7) 四列 multiselect ----------
    labels = {
        "binary":     "BIN_VARS [二元变量]",
        "ordinal":    "ORD_VARS [有序变量]",
        "nominal":    "NOM_VARS [名义变量]",
        "continuous": "CONT_VARS [连续变量]",
    }
    cols = st.columns(4)

    for cat, col in zip(cats, cols):
        with col:
            st.multiselect(
                label      = labels[cat],
                options    = all_vars,                     # 给全集
                key        = f"sel_{cat}",
                on_change  = partial(handle_change, cat),  # 回调里去重
            )
    state.id_col_val   = state.id_col
    state.dep_var_val  = state.dep_var
    state.binary_dv_val = state.binary_dv    

    # ---------- 8) 下一步 ----------
    st.button("下一步: CE1 抽样", on_click=go_to, args=(3,))
# -------- Step 3: CE1 抽样 --------
elif state.step == 3:
    st.header("Step 3: CE1 抽样", divider=True)
    df = state.get("df")
    if df is None:
        st.warning("请先完成前面步骤")
    else:
        assigned   = state.get("assigned", {})
        cont_vars  = assigned.get("continuous", [])
        nom_vars   = assigned.get("nominal", [])
        bin_vars   = assigned.get("binary", [])
        ord_vars   = assigned.get("ordinal", [])
        id_col     = state.get("id_col_val")
        dep_var    = state.get("dep_var_val")
        binary_dv  = state.get("binary_dv_val")
        left, right = st.columns([1, 2], border=True)
        with left:
            with st.form("ce1_form", border=False):
                st.subheader("General & Sampling Variables", divider=True)
                
                label_out_dir = '''**PATH_OUTPUT**  
                :blue-background[Your output destination folder]'''
                out_dir = st.text_input(label_out_dir, os.path.join(BASE_DIR, "streamlit_output"))
                
                label_split_portion = '''**SPLIT_PORTION**  
                :blue-background[Your split portion for modeling portion]'''
                split_portion = st.number_input(label_split_portion, value=0.7)
                
                label_exclusion_if = '''**EXCLUTION_IF**  
                :blue-background[Data exclusion condition; set an actual condition if needed, e.g., inds[dep_var].isna()]'''
                exclusion_if = st.text_input(label_exclusion_if, value="inds[dep_var].isna()")
                
                with st.expander("Optional Macro Variables: These all have defaults that can be used", expanded=False):

                    label_bootstrap = '''**BOOTSTRAP**  
                    :blue-background[Request oversampling, bootstrap if responders size is small (Y/N). Binary DV only]'''
                    bootstrap = st.selectbox(label_bootstrap, ["Y", "N"], index=0)
                    
                    label_oversampled_rr = '''**OVERSAMPLED_RR**  
                    :blue-background[Your desired oversample response rate on the modeling part of data, 0.1-0.5]'''
                    oversampled_rr = st.number_input(label_oversampled_rr, value=0.1)
                    
                    label_min_num_resp = '''**MIN_NUM_RESP**  
                    :blue-background[Minimum number of responders if sampling]'''
                    min_num_resp = st.number_input(label_min_num_resp, value=2000)

                    label_seed = '''**SEED**  
                    :blue-background[Random selection seed if sampling]'''
                    seed = st.number_input(label_seed, value=123456)

                submitted = st.form_submit_button("Run CE1 Sampling")
            if submitted:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log1 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
                lst1 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")
                ce1_config = {
                    "cont_vars": cont_vars,
                    "nom_vars": nom_vars,
                    "bin_vars": bin_vars,
                    "ord_vars": ord_vars,
                    "path_output": out_dir,
                    "inds": df,
                    "id": id_col,
                    "dep_var": dep_var,
                    "binary_dv": binary_dv,
                    "weight": None,
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
                state.ce1_config = ce1_config
                with printto(log=log1, lst=lst1) as logger:
                    sampled_df, rate = CE_Sampling(config=ce1_config, logger=logger)
                state.sampled_df = sampled_df
                state.ce1_rate = rate
                state.log1 = log1
                state.lst1 = lst1
        with right:
            if "sampled_df" in state:
                st.subheader("抽样结果", divider=True)
                st.write(f"样本率：{state.ce1_rate}")
                st.dataframe(state.sampled_df.head(100))
                st.download_button(
                    "下载sampled_df",
                    state.sampled_df.to_csv(index=False),
                    "CE2_sampled_df.csv",
                )
                with st.expander("查看 CE1 日志"):
                    st.code(open(state.log1, encoding="utf-8").read())
                st.button("下一步: CE2 重编码", on_click=go_to, args=(4,))
            else:
                st.info("请在左侧填写参数并运行 CE1。")

# -------- Step 4: CE2 重编码 --------
elif state.step == 4:
    st.header("Step 4: CE2 重编码", divider=True)
    sampled = state.get("sampled_df")
    if sampled is None:
        st.warning("请先完成 CE1 抽样")
    else:
        left, right = st.columns([1, 2], border=True)
        with left:
            with st.form("ce2_form", border=False):
                st.subheader("General & Sampling Variables", divider=True)
                
                label_out_dir = '''**PATH_OUTPUT**  
                :blue-background[Your output destination folder]'''
                out_dir = st.text_input(label_out_dir, os.path.join(BASE_DIR, "streamlit_output"))

                with st.expander("Options are ER for Equal Response method", expanded=False):
                    label_profiling = '''**PROFILING**  
                    :blue-background[Request profiling report on all variables (Y/N)]'''
                    profiling = st.pills(label_profiling, ["Y", "N"], default="Y", key='profiling')

                    label_missrate = '''**MISSRATE**  
                    :blue-background[Maximum missing rate allowed]'''
                    missrate = st.number_input(label_missrate, value=0.75)

                    label_concrate = '''**CONCRATE**  
                    :blue-background[Maximum amount of file that can be in a single value of an independent variable]'''
                    concrate = st.number_input(label_concrate, value=0.9)
                    
                    label_valcnt = '''**VALCNT**  
                    :blue-background[Maximum number of unique values allowed. Set to 0 to allow any number]'''
                    valcnt = st.number_input(label_valcnt, value=50)

                    label_minbinnc = '''**MINBINNC**  
                    :blue-background[Minimum count in a bin to be usable. Set to 0 to use minbinpct only]'''
                    minbinnc = st.number_input(label_minbinnc, value=500)

                    label_minbinnp = '''**MINBINNP**  
                    :blue-background[Minimum percent of file in a bin to be usable. Set to 0 to use minbincnt only]'''
                    minbinnp = st.number_input(label_minbinnp, value=0.05)

                    label_talpha = '''**TALPHA**  
                    :blue-background[T-Test significance level for collapse of bins]'''
                    talpha = st.number_input(label_talpha, value=0.05)

                    label_bonfer = '''**BONFER**  
                    :blue-background[Do Bonferroni adjustment for talpha? (Y/N)]'''
                    bonfer = st.pills(label_bonfer, ["Y", "N"], default="N", key='bonfer')

                    label_nom_method = '''**NOM_METHOD**  
                    :blue-background[Recoding method for nominal variables: Binary, Index or Mean]'''
                    nom_method = st.text_input(label_nom_method, "INDEX")

                    label_pvalue = '''**PVALUE**  
                    :blue-background[P-value threshold to include variable]'''
                    pvalue = st.number_input(label_pvalue, value=0.05)

                    label_min_size = '''**MIN_SIZE**  
                    :blue-background[Minimum missing group size for Equal Response imputation]'''
                    min_size = st.number_input(label_min_size, value=500)

                    label_num_category = '''**NUM_CATEGORY**  
                    :blue-background[Maximum number of categories for profiling variables]'''
                    num_category = st.number_input(label_num_category, value=10)

                    label_equal_dist = '''**EQUAL_DIST**  
                    :blue-background[Use equal distance for dividing variables into groups for profiling (Y/N)]'''
                    equal_dist = st.pills(label_equal_dist, ["Y", "N"], default="N", key='equal_dist')

                    label_p_lo = '''**P_LO**  
                    :blue-background[Lower percentile for checking constant value]'''
                    p_lo = st.number_input(label_p_lo, value=1)

                    label_p_hi = '''**P_HI**  
                    :blue-background[Upper percentile for checking constant value]'''
                    p_hi = st.number_input(label_p_hi, value=99)

                    label_impmethodC = '''**IMPMETHODC**  
                    :blue-background[What method to use for missing imputation for continuous variables?]'''
                    impmethodC = st.text_input(label_impmethodC, "median")

                    label_stdmethodC = '''**STDMETHODC**  
                    :blue-background[Standardization options: any method allowed in proc stdize or NO to skip]'''
                    stdmethodC = st.text_input(label_stdmethodC, "STD")

                    label_cap_flrC = '''**CAP_FLRC**  
                    :blue-background[What method to use for missing imputation for continuous variables?]'''
                    cap_flrC = st.pills(label_cap_flrC, ["Y", "N"], default="Y", key='cap_flrC')

                    label_transformationC = '''**TRANSFORMATIONC**  
                    :blue-background[Include transformed variables in evaluation (Y/N)]'''
                    transformationC = st.pills(label_transformationC, ["Y", "N"], default="Y", key='')

                    label_impmethodO = '''**IMPMETHODO**  
                    :blue-background[What method to use for missing imputation for ordinal variables?]'''
                    impmethodO = st.text_input(label_impmethodO, "mean")

                    label_stdmethodO = '''**STDMETHODO**  
                    :blue-background[Standardization options: any method allowed in proc stdize or NO to skip]'''
                    stdmethodO = st.text_input(label_stdmethodO, "No")

                    label_cap_flrO = '''**CAP_FLRO**  
                    :blue-background[Do you want to do capping/flooring to handle outliers? (Y/N)]'''
                    cap_flrO = st.pills(label_cap_flrO, ["Y", "N"], default="N", key='cap_flrO')

                    label_transformationO = '''**TRANSFORMATIONO**  
                    :blue-background[Include transformed variables in evaluation? (Y/N)]'''
                    transformationO = st.pills(label_transformationO, ["Y", "N"], default="N", key='transformationO')

                submitted2 = st.form_submit_button("运行 CE2")
            if submitted2:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log2 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
                lst2 = os.path.join(log_dir, "02_CE_EDA_REcode_LST_File.lst")
                ce2_cfg = state.ce1_config.copy()
                ce2_cfg.update(
                    {
                        "inds": sampled,
                        "profiling": profiling,
                        "pvalue": pvalue,
                        "min_size": min_size,
                        "num_category": num_category,
                        "equal_dist": equal_dist,
                    }
                )
                with printto(log=log2, lst=lst2) as logger:
                    rec, prof = CE_EDA_Recode(
                        indsn=sampled, config=ce2_cfg, logger=logger
                    )
                state.recoded_df = rec
                state.profile_df = prof
                state.log2 = log2
                state.lst2 = lst2
        with right:
            if "recoded_df" in state:
                st.subheader("重编码结果", divider=True)
                st.dataframe(state.recoded_df.head(100))
                st.download_button(
                    "下载重编码数据",
                    state.recoded_df.to_csv(index=False),
                    "CE2_Recoded.csv",
                )
                st.subheader("Profiling", divider=True)
                st.dataframe(state.profile_df)
                st.download_button(
                    "下载 Profiling 报表",
                    state.profile_df.to_csv(index=False),
                    "Profile.csv",
                )
                with st.expander("查看 CE2 日志"):
                    st.code(open(state.log2, encoding="utf-8").read())
            else:
                st.info("请在左侧填写参数并运行 CE2。")
