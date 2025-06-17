# -*- coding: utf-8 -*-
# Auto-generated ordinal recode

# --- Recode variable: Pclass ---
df['R1_Pclass'] = df['Pclass'].fillna(2.2667731629392973)
df['R1_Pclass'].attrs['label'] = 'Pclass: Recode (-)'

# --- Final clip & label for Pclass ---
df['R1_Pclass'] = df['R1_Pclass'].clip(lower=1.0, upper=3.0)
df['R1_Pclass'].attrs['label'] = 'Pclass: Recode (-)'
df['Pclass'] = df['R1_Pclass']
df['Pclass'].attrs['label'] = 'Pclass Pclass (-)'


# KEEP_LIST_O
KEEP_LIST_O = ['Pclass']
