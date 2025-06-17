# -*- coding: utf-8 -*-
# Auto-generated nominal recode

# Nominal recode for Embarked
df['R1_Embarked'] = df['Embarked'].map({'S': 0.32397408207343414, 'Q': 0.46808510638297873, 'C': 0.5877192982456141, nan: nan}).fillna(0)

# KEEP_LIST_N
KEEP_LIST_N = ['R1_Embarked']
