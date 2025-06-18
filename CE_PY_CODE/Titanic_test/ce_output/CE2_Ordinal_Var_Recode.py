# -*- coding: utf-8 -*-
# Auto-generated ordinal recode

# --- Recode variable: Pclass ---
df['R1_Pclass'] = df['Pclass'].fillna(2.47695)
df['R1_Pclass'].attrs['label'] = 'Pclass: Recode (-)'

# --- Final clip & label for Pclass ---
df['R1_Pclass'] = df['R1_Pclass'].clip(lower=1.0, upper=3.0)
df['R1_Pclass'].attrs['label'] = 'Pclass: Recode (-)'
df['Pclass'] = df['R1_Pclass']
df['Pclass'].attrs['label'] = 'Pclass Pclass (-)'

# --- Recode variable: Parch ---
df['R1_Parch'] = df['Parch'].fillna(0.3293)
df['R1_Parch'].attrs['label'] = 'Parch: Recode (+)'

# --- Final clip & label for Parch ---
df['R1_Parch'] = df['R1_Parch'].clip(lower=0.0, upper=6.0)
df['R1_Parch'].attrs['label'] = 'Parch: Recode (+)'
df['Parch'] = df['R1_Parch']
df['Parch'].attrs['label'] = 'Parch Parch (+)'


# KEEP_LIST_O
KEEP_LIST_O = ['Pclass', 'Parch']
