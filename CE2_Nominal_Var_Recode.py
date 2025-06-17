# -*- coding: utf-8 -*-
# Auto-generated nominal recode

# Nominal recode for Embarked
df['R1_Embarked_X1'] = df['Embarked'].isin(['S', 'Q']).astype(int)

# KEEP_LIST_N
KEEP_LIST_N = ['R1_Embarked_X1']
