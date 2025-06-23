import os
import sys
import pandas as pd
import numpy as np
import streamlit as st

# 将 CE_PY_CODE 目录加入搜索路径
BASE_DIR = os.path.dirname(__file__)
sys.path.append(os.path.join(BASE_DIR, "CE_PY_CODE"))

from CE1_Sampling import CE_Sampling
from CE2_EDA_Recode import CE_EDA_Recode
from CE_Log_Tool import printto

# Streamlit 页面配置
st.set_page_config(page_title="CeFlow", layout="wide")

# 初始化 session_state
state = st.session_state
if 'step' not in state:
    state.step = 1

# 侧边栏导航 标签
step_labels = ["1. 数据上传", "2. 变量确认", "3. CE1 抽样", "4. CE2 重编码"]
# 渲染 radio，key 与 session_state 绑定，无需 index
sel = st.sidebar.radio("选择步骤", step_labels, key='sidebar_sel')
# 根据侧边栏选择更新 step
state.step = step_labels.index(sel) + 1

# 回调：按钮点击更新 sidebar_sel

def go_to(step):
    state.sidebar_sel = step_labels[step-1]

st.title("CeFlow Streamlit IDE")

# -------- Step 1: 上传和预览数据 --------
if state.step == 1:
    st.header("Step 1: 上传和预览数据")
    uploaded = st.file_uploader("上传数据 (csv/xlsx)", type=["csv", "xlsx", "xls"])
    path_input = st.text_input("或输入文件路径")
    df = None
    if uploaded:
        df = pd.read_csv(uploaded) if uploaded.name.lower().endswith('.csv') else pd.read_excel(uploaded)
    elif path_input:
        df = pd.read_csv(path_input) if path_input.lower().endswith('.csv') else pd.read_excel(path_input)
    if df is not None:
        state.df = df
        st.subheader("数据预览 (仅前 100 行)")
        st.dataframe(df.head(100))
        st.button("下一步: 变量确认", on_click=go_to, args=(2,))

# -------- Step 2: 变量确认与配置 --------
elif state.step == 2:
    st.header("Step 2: 变量确认与配置")
    df = state.get('df')
    if df is None:
        st.warning("请先完成 Step 1 上传数据")
    else:
        st.subheader("数据预览 (前 5 行)")
        st.dataframe(df.head())
        st.subheader("标识 & 因变量配置")
        state.id_col    = st.selectbox("ID 列", df.columns.tolist(), index=0)
        state.dep_var   = st.selectbox("因变量", df.columns.tolist(), index=1)
        state.binary_dv = st.selectbox("因变量类型 (Binary DV)", ["Y","N"], index=0)

        # 自动检测变量类型
        num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        auto_bin  = [c for c in num_cols if df[c].nunique() == 2]
        auto_cont = [c for c in num_cols if c not in auto_bin]
        obj_cols  = df.select_dtypes(include=['object','category']).columns.tolist()
        auto_nom  = obj_cols.copy()
        auto_ord  = []

        st.subheader("变量类型分配 (可手动调整)")
        default_bin  = state.get('bin_vars', auto_bin)
        default_nom  = state.get('nom_vars', auto_nom)
        default_ord  = state.get('ord_vars', auto_ord)
        default_cont = state.get('cont_vars', auto_cont)

        state.bin_vars  = st.multiselect("二元变量", df.columns.tolist(), default=default_bin)
        state.nom_vars  = st.multiselect("名义变量", df.columns.tolist(), default=default_nom)
        state.ord_vars  = st.multiselect("有序变量", df.columns.tolist(), default=default_ord)
        state.cont_vars = st.multiselect("连续变量", df.columns.tolist(), default=default_cont)

        st.button("下一步: CE1 抽样", on_click=go_to, args=(3,))

# -------- Step 3: CE1 抽样 --------
elif state.step == 3:
    st.header("Step 3: CE1 抽样")
    df = state.get('df')
    if df is None:
        st.warning("请先完成前面步骤")
    else:
        left, right = st.columns([1, 2])
        with left:
            with st.form("ce1_form"):
                st.subheader("抽样参数配置")
                split_portion   = st.number_input("建模样本比例", value=0.7)
                bootstrap       = st.selectbox("Bootstrap", ["Y","N"], index=0)
                oversampled_rr  = st.number_input("目标响应率", value=0.1)
                min_num_resp    = st.number_input("最小响应数", value=2000)
                seed            = st.number_input("随机种子", value=123456)
                out_dir         = st.text_input("输出目录", os.path.join(BASE_DIR, "streamlit_output"))
                submitted       = st.form_submit_button("运行 CE1")
            if submitted:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log1 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
                lst1 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")
                ce1_config = {"cont_vars": state.cont_vars, "nom_vars": state.nom_vars,
                              "bin_vars": state.bin_vars, "ord_vars": state.ord_vars,
                              "path_output": out_dir, "inds": df, "id": state.id_col,
                              "dep_var": state.dep_var, "binary_dv": state.binary_dv,
                              "weight": None, "split_portion": split_portion,
                              "exclusion_if": "inds[dep_var].isna()", "DS_present": "N",
                              "prefix": "R1_", "keep_list": [], "bootstrap": bootstrap,
                              "oversampled_rr": oversampled_rr, "min_num_resp": min_num_resp,
                              "seed": int(seed)}
                state.ce1_config = ce1_config
                with printto(log=log1, lst=lst1) as logger:
                    sampled_df, rate = CE_Sampling(config=ce1_config, logger=logger)
                state.sampled_df = sampled_df
                state.ce1_rate  = rate
                state.log1      = log1
                state.lst1      = lst1
        with right:
            if 'sampled_df' in state:
                st.subheader("抽样结果")
                st.write(f"样本率：{state.ce1_rate}")
                st.dataframe(state.sampled_df)
                with st.expander("查看 CE1 日志"):
                    st.code(open(state.log1, encoding="utf-8").read())
                st.button("下一步: CE2 重编码", on_click=go_to, args=(4,))
            else:
                st.info("请在左侧填写参数并运行 CE1。")

# -------- Step 4: CE2 重编码 --------
elif state.step == 4:
    st.header("Step 4: CE2 重编码")
    sampled = state.get('sampled_df')
    if sampled is None:
        st.warning("请先完成 CE1 抽样")
    else:
        left, right = st.columns([1, 2])
        with left:
            with st.form("ce2_form"):
                profiling    = st.selectbox("生成 Profiling 报告", ["Y","N"], index=0)
                pvalue       = st.number_input("显著性水平", value=0.05)
                min_size     = st.number_input("最小缺失组大小", value=500)
                num_category = st.number_input("分箱数", value=10)
                equal_dist   = st.selectbox("等距分箱", ["N","Y"], index=0)
                out_dir      = st.text_input("输出目录", os.path.join(BASE_DIR, "streamlit_output"))
                submitted2   = st.form_submit_button("运行 CE2")
            if submitted2:
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log2 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
                lst2 = os.path.join(log_dir, "02_CE_EDA_REcode_LST_File.lst")
                ce2_cfg = state.ce1_config.copy()
                ce2_cfg.update({"inds": sampled, "profiling": profiling,
                                "pvalue": pvalue, "min_size": min_size,
                                "num_category": num_category, "equal_dist": equal_dist})
                with printto(log=log2, lst=lst2) as logger:
                    rec, prof = CE_EDA_Recode(indsn=sampled, config=ce2_cfg, logger=logger)
                state.recoded_df  = rec
                state.profile_df  = prof
                state.log2        = log2
                state.lst2        = lst2
        with right:
            if 'recoded_df' in state:
                st.subheader("重编码结果")
                st.dataframe(state.recoded_df)
                st.download_button("下载重编码数据", state.recoded_df.to_csv(index=False), "CE2_Recoded.csv")
                st.dataframe(state.profile_df)
                st.download_button("下载 Profiling 报表", state.profile_df.to_csv(index=False), "Profile.csv")
                with st.expander("查看 CE2 日志"):
                    st.code(open(state.log2, encoding="utf-8").read())
            else:
                st.info("请在左侧填写参数并运行 CE2。")
