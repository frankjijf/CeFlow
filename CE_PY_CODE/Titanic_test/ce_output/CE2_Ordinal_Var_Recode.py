# -*- coding: utf-8 -*-
# Auto-generated ordinal recode

# --- Recode variable: Pclass ---
df.loc[df['Pclass'].isna(), 'R1_Pclass'] = 3
df.loc[df['Pclass'].notna(), 'R1_Pclass'] = df['Pclass']

# --- Recode variable: Parch ---
df.loc[df['Parch'].isna(), 'R1_Parch'] = 0
df.loc[df['Parch'].notna(), 'R1_Parch'] = df['Parch']


# KEEP_LIST_O
KEEP_LIST_O = ['R1_Pclass', 'R1_Parch']
