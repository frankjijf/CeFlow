# -*- coding: utf-8 -*-
# Auto-generated nominal recode

# Recode variable: Sex
df['R1_Sex'] = df['Sex'].map({'male': 40.18883160444748, 'female': 346.83414509100226}).fillna(0)
# Recode variable: Embarked
df['R1_Embarked'] = df['Embarked'].map({'S': 77.83548118227864, 'Q': 134.53536754507627, 'C': 199.93245525160418, nan: 199.93245525160418}).fillna(0)

# KEEP_LIST_N
KEEP_LIST_N = ['R1_Sex', 'R1_Embarked']
