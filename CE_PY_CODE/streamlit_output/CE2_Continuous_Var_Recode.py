# -*- coding: utf-8 -*-
# Auto-generated continuous recode

# --- Recode variable: Age ---
df.loc[df['Age'].isna(), 'R1_Age'] = 28.0
df.loc[df['Age'].notna(), 'R1_Age'] = df['Age']
df['R1_LN_Age'] = np.log(np.maximum(df['R1_Age'], 1e-5))

# --- Recode variable: Fare ---
df.loc[df['Fare'].isna(), 'R1_Fare'] = 11.5
df.loc[df['Fare'].notna(), 'R1_Fare'] = df['Fare']
df['R1_SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))


# KEEP_LIST_C
KEEP_LIST_C = ['R1_Age', 'R1_Fare']
