# -*- coding: utf-8 -*-
# Auto-generated continuous recode

# --- Recode variable: PassengerId ---
df.loc[df['PassengerId'].isna(), 'R1_PassengerId'] = 435
df.loc[df['PassengerId'].notna(), 'R1_PassengerId'] = df['PassengerId']
df['R1_LN_PassengerId'] = np.log(np.maximum(df['R1_PassengerId'], 1e-5))

# --- Recode variable: Pclass ---
df.loc[df['Pclass'].isna(), 'R1_Pclass'] = 3
df.loc[df['Pclass'].notna(), 'R1_Pclass'] = df['Pclass']

# --- Recode variable: Age ---
df.loc[df['Age'].isna(), 'R1_Age'] = 28.0
df.loc[df['Age'].notna(), 'R1_Age'] = df['Age']
df['R1_LN_Age'] = np.log(np.maximum(df['R1_Age'], 1e-5))

# --- Recode variable: SibSp ---
df.loc[df['SibSp'].isna(), 'R1_SibSp'] = 0
df.loc[df['SibSp'].notna(), 'R1_SibSp'] = df['SibSp']
df['R1_LN_SibSp'] = np.log(np.maximum(df['R1_SibSp'], 1e-5))

# --- Recode variable: Parch ---
df.loc[df['Parch'].isna(), 'R1_Parch'] = 0
df.loc[df['Parch'].notna(), 'R1_Parch'] = df['Parch']
df['R1_LN_Parch'] = np.log(np.maximum(df['R1_Parch'], 1e-5))

# --- Recode variable: Fare ---
df.loc[df['Fare'].isna(), 'R1_Fare'] = 11.5
df.loc[df['Fare'].notna(), 'R1_Fare'] = df['Fare']
df['R1_SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))


# KEEP_LIST_C
KEEP_LIST_C = ['R1_PassengerId', 'R1_Pclass', 'R1_Age', 'R1_SibSp', 'R1_Parch', 'R1_Fare']
