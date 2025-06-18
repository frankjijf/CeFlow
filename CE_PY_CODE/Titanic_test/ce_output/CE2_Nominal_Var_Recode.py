# -*- coding: utf-8 -*-
# Auto-generated nominal recode

# Nominal recode for Embarked
df['R1_Embarked'] = df['Embarked'].map({'S': 0.07783548118227865, 'Q': 0.13453536754507628, 'C': 0.1953125, nan: nan}).fillna(0)
# Nominal recode for Cabin2
df['R1_Cabin2'] = df['Cabin2'].map({'E': 0.0, 'S': 0.07629281396910678, nan: nan, 'Q': 0.1295577967416602, 'B': 0.1509433962264151, 'C': 0.19806403574087864, 'A': 1.0, 'D': 1.0, 'G': 1.0}).fillna(0)

# KEEP_LIST_N
KEEP_LIST_N = ['R1_Embarked', 'R1_Cabin2']
