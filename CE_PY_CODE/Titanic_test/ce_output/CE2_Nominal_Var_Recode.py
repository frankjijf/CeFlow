# -*- coding: utf-8 -*-
# Auto-generated nominal recode

# Recode variable: Embarked
df['R1_Embarked'] = df['Embarked'].map({'S': 77.83548118227864, 'Q': 134.53536754507627, 'C': 199.93245525160418, nan: 199.93245525160418}).fillna(0)

# KEEP_LIST_N
KEEP_LIST_N = ['R1_Embarked']
