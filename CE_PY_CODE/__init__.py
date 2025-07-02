"""Core CeFlow modules with lazy imports."""

__all__ = ["CE_Sampling", "CE_EDA_Recode", "CE_Var_Redu", "printto"]


def __getattr__(name):
    if name == "CE_Sampling":
        from .ce1_sampling import CE_Sampling

        return CE_Sampling
    if name == "CE_EDA_Recode":
        from .ce2_eda_recode import CE_EDA_Recode

        return CE_EDA_Recode
    if name == "CE_Var_Redu":
        from .ce3_var_redu import CE_Var_Redu

        return CE_Var_Redu
    if name == "printto":
        from .ce_log_tool import printto

        return printto
    raise AttributeError(f"module {__name__} has no attribute {name}")
