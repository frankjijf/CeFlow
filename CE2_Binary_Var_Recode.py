# -*- coding: utf-8 -*-
# Auto-generated binary recode

# Recode variable: IsFemale
df['R1_IsFemale'] = (df['IsFemale']==1).astype(int)
df['R1_IsFemale'].attrs['label'] = 'Is Female: Binary Recode'

# Recode variable: IsChild
df['R1_IsChild'] = (df['IsChild']==1).astype(int)
df['R1_IsChild'].attrs['label'] = 'Is Child: Binary Recode'


# KEEP_LIST_B
KEEP_LIST_B = ['R1_IsFemale', 'R1_IsChild']
