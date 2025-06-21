# -*- coding: utf-8 -*-
# Auto-generated ordinal recode

# --- Recode variable: Pclass ---
df.loc[df['Pclass'].isna(), 'R1_Pclass'] = 3
df.loc[df['Pclass'].notna(), 'R1_Pclass'] = df['Pclass']
df['R1_Pclass'] = df['R1_Pclass'].clip(lower=1, upper=3)

# --- Recode variable: Parch ---
df.loc[df['Parch'].isna(), 'R1_Parch'] = 0
df.loc[df['Parch'].notna(), 'R1_Parch'] = df['Parch']
df['R1_Parch'] = df['R1_Parch'].clip(lower=0, upper=6.0)
df['R1_LN_Parch'] = np.log(np.maximum(df['R1_Parch'], 1e-5))


# KEEP_LIST_O
KEEP_LIST_O = ['R1_Pclass', 'R1_LN_Parch']
