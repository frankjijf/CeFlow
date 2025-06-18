# -*- coding: utf-8 -*-
# Auto-generated continuous recode

# --- Recode variable: Age ---
df['R1_Age'] = df['Age'].fillna(28.0)
df['R1_Age'] = df['R1_Age'].clip(lower=0.42, upper=71.0)
df['R1_Age'].attrs['label'] = 'Age: Recode (-)'
df['IV_Age'] = -1/np.maximum(df['R1_Age'], 1e-5)
df['IV_Age'].attrs['label'] = 'Age IV_ (-)'

# --- Recode continuous variable: Age ---
df['R1_Age'] = df['Age'].fillna(28.0)
df['R1_Age'] = df['R1_Age'].clip(lower=0.42, upper=71.0)
df['R1_Age'].attrs['label'] = 'Age: Recode (-)'
df['IV_Age'] = -1/np.maximum(df['R1_Age'], 1e-5)
df['IV_Age'].attrs['label'] = 'Age IV (-)'

# --- Recode variable: Fare ---
df['R1_Fare'] = df['Fare'].fillna(11.5)
df['R1_Fare'] = df['R1_Fare'].clip(lower=0.0, upper=234.02505000000002)
df['R1_Fare'].attrs['label'] = 'Fare: Recode (+)'
df['SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))
df['SR_Fare'].attrs['label'] = 'Fare SR_ (+)'

# --- Recode continuous variable: Fare ---
df['R1_Fare'] = df['Fare'].fillna(11.5)
df['R1_Fare'] = df['R1_Fare'].clip(lower=0.0, upper=234.02505000000002)
df['R1_Fare'].attrs['label'] = 'Fare: Recode (+)'
df['SR_Fare'] = np.sqrt(np.maximum(df['R1_Fare'], 0))
df['SR_Fare'].attrs['label'] = 'Fare SR (+)'


# KEEP_LIST_C
KEEP_LIST_C = ['IV_Age', 'SR_Fare']
