# -*- coding: utf-8 -*-
# Auto-generated continuous recode

# --- Recode variable: Age ---
df['R1_Age'] = df['Age'].fillna(28.0)
df['R1_Age'] = df['R1_Age'].clip(lower=0.42, upper=79.805)
df['R1_Age'].attrs['label'] = 'Age: Recode (-)'
df['LN_Age'] = np.log(np.maximum(df['R1_Age'], 1e-5))
df['LN_Age'].attrs['label'] = 'Age LN_ (-)'

# --- Recode continuous variable: Age ---
df['R1_Age'] = df['Age'].fillna(28.0)
df['R1_Age'] = df['R1_Age'].clip(lower=0.42, upper=79.805)
df['R1_Age'].attrs['label'] = 'Age: Recode (-)'
df['LN_Age'] = np.log(np.maximum(df['R1_Age'], 1e-5))
df['LN_Age'].attrs['label'] = 'Age LN (-)'

# --- Recode variable: Fare ---
df['R1_Fare'] = df['Fare'].fillna(14.4542)
df['R1_Fare'] = df['R1_Fare'].clip(lower=0.0, upper=358.00933000000055)
df['R1_Fare'].attrs['label'] = 'Fare: Recode (+)'
df['SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))
df['SR_Fare'].attrs['label'] = 'Fare SR_ (+)'

# --- Recode continuous variable: Fare ---
df['R1_Fare'] = df['Fare'].fillna(14.4542)
df['R1_Fare'] = df['R1_Fare'].clip(lower=0.0, upper=358.00933000000055)
df['R1_Fare'].attrs['label'] = 'Fare: Recode (+)'
df['SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))
df['SR_Fare'].attrs['label'] = 'Fare SR (+)'


# KEEP_LIST_C
KEEP_LIST_C = ['LN_Age', 'SR_Fare']
