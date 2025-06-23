# -*- coding: utf-8 -*-
# Auto-generated ordinal recode

# --- Recode variable: Pclass ---
df.loc[df['Pclass'].isna(), 'R1_Pclass'] = 3
df.loc[df['Pclass'].notna(), 'R1_Pclass'] = df['Pclass']


# KEEP_LIST_O
KEEP_LIST_O = ['R1_Pclass']
