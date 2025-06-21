# -*- coding: utf-8 -*-
# Auto-generated binary recode

# Recode variable: IsFemale
df['R1_IsFemale'] = (df['IsFemale'] == 1).fillna(False).astype(int)
# Recode variable: IsChild
df['R1_IsChild'] = (df['IsChild'] == 1).fillna(False).astype(int)

# KEEP_LIST_B
KEEP_LIST_B = ['R1_IsFemale', 'R1_IsChild']
