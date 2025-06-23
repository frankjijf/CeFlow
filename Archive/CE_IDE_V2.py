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
if 'ce1_done' not in state:
    state.ce1_done = False
if 'ce2_done' not in state:
    state.ce2_done = False

if 'step' not in st.session_state:
    st.session_state.step = 1

# —— 下一步回调函数 —— 
def go_next(step):
    st.session_state.step = step

# 侧边栏步骤导航
step_labels = ["1. 数据上传", "2. 变量确认", "3. CE1 抽样", "4. CE2 重编码"]
sel = st.sidebar.radio("选择步骤", step_labels, index=state.step-1)
state.step = step_labels.index(sel) + 1

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
        st.write("数据预览", df.head())
        state.df = df
        if st.button("下一步: 变量确认"):
            state.step = 2

# -------- Step 2: 自动判断变量属性 --------
elif state.step == 2:
    st.header("Step 2: 变量确认与配置")
    df = state.get('df')
    if df is None:
        st.warning("请先完成 Step 1 上传数据")
    else:
        # 展示预览
        st.write("数据预览", df.head(100))

        # ID、因变量及类型
        st.subheader("标识 & 因变量配置")
        state.id_col = st.selectbox("ID 列", options=df.columns.tolist(), index=0)
        state.dep_var = st.selectbox("因变量", options=df.columns.tolist(), index=1)
        state.binary_dv = st.selectbox("因变量类型 (Binary DV)", ["Y", "N"], index=0)

        # 自动检测变量类型
        num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        auto_bin = [c for c in num_cols if df[c].nunique() == 2]
        auto_cont = [c for c in num_cols if c not in auto_bin]
        obj_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()
        auto_nom = [c for c in obj_cols if c not in auto_bin]
        auto_ord = []

        st.subheader("变量类型分配 (可手动调整)")
        state.bin_vars = st.multiselect("二元变量", options=df.columns.tolist(), default=auto_bin)
        state.nom_vars = st.multiselect("名义变量", options=df.columns.tolist(), default=auto_nom)
        state.ord_vars = st.multiselect("有序变量", options=df.columns.tolist(), default=auto_ord)
        state.cont_vars = st.multiselect("连续变量", options=df.columns.tolist(), default=auto_cont)

        if st.button("下一步: CE1 抽样"):
            state.step = 3

# -------- Step 3: CE1 抽样 --------
elif state.step == 3:
    st.header("Step 3: CE1 抽样")
    df = state.get('df')
    if df is None:
        st.warning("请先完成前面步骤")
    elif not state.ce1_done:
        left, right = st.columns([1, 2])
        with left:
            st.subheader("抽样参数配置")
            split_portion = st.number_input("建模样本比例", value=0.7)
            bootstrap = st.selectbox("Bootstrap", ["Y", "N"], index=0)
            oversampled_rr = st.number_input("目标响应率", value=0.1)
            min_num_resp = st.number_input("最小响应数", value=2000)
            seed = st.number_input("随机种子", value=123456)
            out_dir = st.text_input("输出目录", os.path.join(BASE_DIR, "streamlit_output"))

            if st.button("运行 CE1"):
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log1 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
                lst1 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")

                # 构造 CE1 配置
                ce1_config = {
                    "cont_vars": state.cont_vars,
                    "nom_vars": state.nom_vars,
                    "bin_vars": state.bin_vars,
                    "ord_vars": state.ord_vars,
                    "path_output": out_dir,
                    "inds": df,
                    "id": state.id_col,
                    "dep_var": state.dep_var,
                    "binary_dv": state.binary_dv,
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
                # 保存 CE1 配置
                state.ce1_config = ce1_config

                with printto(log=log1, lst=lst1) as logger:
                    sampled_df, ce1_rate = CE_Sampling(config=ce1_config, logger=logger)

                state.sampled_df = sampled_df
                state.ce1_rate = ce1_rate
                state.log1 = log1
                state.lst1 = lst1
                state.ce1_done = True

        with right:
            st.subheader("运行提示")
            st.write("在左侧配置抽样参数，点击运行。完成后在此查看结果。")
    else:
        st.subheader("抽样结果")
        st.write(f"样本率： {state.ce1_rate}")
        st.dataframe(state.sampled_df.head(100))
        with st.expander("查看 CE1 日志"):
            st.code(open(state.log1, encoding="utf-8").read())
        if st.button("下一步: CE2 重编码"):
            state.step = 4

# -------- Step 4: CE2 重编码 --------
elif state.step == 4:
    st.header("Step 4: CE2 重编码")
    sampled = state.get('sampled_df')
    if sampled is None:
        st.warning("请先完成 CE1 抽样")
    elif not state.ce2_done:
        left, right = st.columns([1, 2])
        with left:
            st.subheader("重编码参数配置")
            profiling = st.selectbox("生成 Profiling 报告", ["Y", "N"], index=0)
            pvalue = st.number_input("显著性水平", value=0.05)
            min_size = st.number_input("最小缺失组大小", value=500)
            num_category = st.number_input("分箱数", value=10)
            equal_dist = st.selectbox("等距分箱", ["N", "Y"], index=0)
            out_dir = st.text_input("输出目录", os.path.join(BASE_DIR, "streamlit_output"))

            if st.button("运行 CE2"):
                os.makedirs(out_dir, exist_ok=True)
                log_dir = os.path.join(out_dir, "log")
                os.makedirs(log_dir, exist_ok=True)
                log2 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
                lst2 = os.path.join(log_dir, "02_CE_EDA_Recode_LST_File.lst")

                # 使用 CE1 配置并新增 CE2 参数
                ce2_config = state.ce1_config.copy()
                ce2_config.update({
                    "inds": sampled,
                    "profiling": profiling,
                    "pvalue": pvalue,
                    "min_size": min_size,
                    "num_category": num_category,
                    "equal_dist": equal_dist,
                })

                with printto(log=log2, lst=lst2) as logger:
                    recoded_df, profile_df = CE_EDA_Recode(indsn=sampled, config=ce2_config, logger=logger)

                state.recoded_df = recoded_df
                state.profile_df = profile_df
                state.log2 = log2
                state.lst2 = lst2
                state.ce2_done = True

        with right:
            st.subheader("运行提示")
            st.write("在左侧配置重编码参数，点击运行。完成后在此查看结果。")
    else:
        st.subheader("重编码结果")
        st.dataframe(state.recoded_df.head(100))
        st.download_button("下载重编码数据", state.recoded_df.to_csv(index=False), "CE2_Recoded.csv")
        st.dataframe(state.profile_df)
        st.download_button("下载 Profiling 报表", state.profile_df.to_csv(index=False), "Profile.csv")
        with st.expander("查看 CE2 日志"):
            st.code(open(state.log2, encoding="utf-8").read()) 