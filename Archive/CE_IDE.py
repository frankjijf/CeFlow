import os
import sys
import pandas as pd
import streamlit as st

# 将 CE_PY_CODE 目录加入搜索路径
BASE_DIR = os.path.dirname(__file__)
sys.path.append(os.path.join(BASE_DIR, "CE_PY_CODE"))

from CE1_Sampling import CE_Sampling
from CE2_EDA_Recode import CE_EDA_Recode
from CE_Log_Tool import printto

st.set_page_config(page_title="CeFlow", layout="wide")
st.title("CeFlow Streamlit 界面")

# ---------------------- 数据读取 ----------------------
st.sidebar.header("数据读取")
uploaded = st.sidebar.file_uploader("上传数据文件 (csv/xlsx)", type=["csv", "xlsx", "xls"])
path_input = st.sidebar.text_input("或输入文件路径")

df = None
if uploaded is not None:
    if uploaded.name.lower().endswith(".csv"):
        df = pd.read_csv(uploaded)
    else:
        df = pd.read_excel(uploaded)
elif path_input:
    if path_input.lower().endswith(".csv"):
        df = pd.read_csv(path_input)
    elif path_input.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path_input)

if df is not None:
    st.write("数据预览", df.head())
    cols = df.columns.tolist()

    # -------------------- 变量列表输入 --------------------
    st.sidebar.header("变量列表")
    id_col = st.sidebar.selectbox("ID 列", cols)
    dep_var = st.sidebar.selectbox("因变量", cols)
    bin_vars = st.sidebar.text_input("二元变量 (逗号分隔)")
    nom_vars = st.sidebar.text_input("名义变量 (逗号分隔)")
    ord_vars = st.sidebar.text_input("有序变量 (逗号分隔)")
    cont_vars = st.sidebar.text_input("连续变量 (逗号分隔)")

    # -------------------- 参数配置 --------------------
    st.sidebar.header("CE1 抽样参数")
    split_portion = st.sidebar.number_input("建模样本比例", value=0.7)
    bootstrap = st.sidebar.selectbox("Bootstrap", ["Y", "N"], index=0)
    oversampled_rr = st.sidebar.number_input("目标响应率", value=0.1)
    min_num_resp = st.sidebar.number_input("最小响应数", value=2000)
    seed = st.sidebar.number_input("随机种子", value=123456)

    st.sidebar.header("CE2 重编码参数")
    profiling = st.sidebar.selectbox("生成 Profiling 报告", ["Y", "N"], index=0)
    pvalue = st.sidebar.number_input("显著性水平", value=0.05)
    min_size = st.sidebar.number_input("最小缺失组大小", value=500)
    num_category = st.sidebar.number_input("分箱数", value=10)
    equal_dist = st.sidebar.selectbox("等距分箱", ["N", "Y"], index=0)

    out_dir = st.sidebar.text_input("输出目录", os.path.join(BASE_DIR, "streamlit_output"))

    if st.sidebar.button("运行流程"):
        with st.spinner("执行 CE1 抽样..."):
            os.makedirs(out_dir, exist_ok=True)
            log_dir = os.path.join(out_dir, "log")
            os.makedirs(log_dir, exist_ok=True)
            log1 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
            lst1 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")

            config = {
                "cont_vars": [v.strip() for v in cont_vars.split(',') if v.strip()],
                "nom_vars": [v.strip() for v in nom_vars.split(',') if v.strip()],
                "bin_vars": [v.strip() for v in bin_vars.split(',') if v.strip()],
                "ord_vars": [v.strip() for v in ord_vars.split(',') if v.strip()],
                "path_output": out_dir,
                "inds": df,
                "id": id_col,
                "dep_var": dep_var,
                "binary_dv": "Y",
                "weight": None,
                "split_portion": split_portion,
                "exclusion_if": "inds[dep_var].isna()",
                "DS_present": "N",
                "prefix": "R1_",
                "path_DS": "/mnt/projects/shared/pst_qmgisi/Modeling/CE/",
                "bootstrap": bootstrap,
                "oversampled_rr": oversampled_rr,
                "min_num_resp": min_num_resp,
                "seed": seed,
                "profiling": profiling,
                "missrate": 0.75,
                "concrate": 0.9,
                "valcnt": 50,
                "minbinnc": 500,
                "minbinnp": 0.05,
                "talpha": 0.05,
                "bonfer": "N",
                "nom_method": "INDEX",
                "pvalue": pvalue,
                "min_size": min_size,
                "num_category": num_category,
                "equal_dist": equal_dist,
                "p_lo": 1,
                "p_hi": 99,
                "impmethodC": "median",
                "stdmethodC": "STD",
                "cap_flrC": "Y",
                "transformationC": "Y",
                "impmethodO": "mean",
                "stdmethodO": "No",
                "cap_flrO": "N",
                "transformationO": "N",
            }
            with printto(log=log1, lst=lst1) as logger:
                sampled, rate = CE_Sampling(config=config, logger=logger)

        with st.spinner("执行 CE2 重编码..."):
            log2 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
            lst2 = os.path.join(log_dir, "02_CE_EDA_Recode_LST_File.lst")
            with printto(log=log2, lst=lst2) as logger:
                recoded, profile = CE_EDA_Recode(sampled, config=config, logger=logger)

        st.success("流程执行完成")
        st.write("样本率表", rate)
        st.write("重编码结果", recoded.head())

        st.download_button("下载重编码数据", recoded.to_csv(index=False), "CE2_Recoded.csv")
        st.download_button("下载 Profile", profile.to_csv(index=False), "Profile.csv")

        log_txt1 = open(log1, encoding="utf-8").read()
        log_txt2 = open(log2, encoding="utf-8").read()
        with st.expander("查看 CE1 日志"):
            st.code(log_txt1)
        with st.expander("查看 CE2 日志"):
            st.code(log_txt2)