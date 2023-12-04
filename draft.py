import pandas as pd

transitions = pd.read_csv('UltimateSplash_Pos_v2.csv')
transition_subset = transitions[transitions['PrecursorName'] == '14:0-13:0-14:0 TG-d5']

transition_subset.to_csv('transition_subset.csv', index=False)

# Add column in transitions dataframe called "Precursor", which contains True and False values. True if the value in ProductName is precursor, False if not
transitions['Precursor'] = transitions['ProductName'].str.contains('precursor')

transitions['Precursor'].value_counts()
# Number of NAs in transitions['Precursor']
transitions['Precursor'].isna().sum()

# Merge 





transition_subset = pd.read_csv('transition_subset.csv')
# Subset transition_subset for columns "PrecursorAdduct", "ProductName", "ProductMz"
transition_subset = transition_subset[['PrecursorName', 'PrecursorAdduct', 'ProductName', 'ProductMz']]
# Remove "1+" from the values in the "PrecursorAdduct" column
transition_subset['PrecursorAdduct'] = transition_subset['PrecursorAdduct'].str.replace('1+', '')
# Remove rows where the value in the "ProductMz" column is NaN
transition_subset = transition_subset.dropna(subset=['ProductMz'])

df_lip_example = pd.read_csv('C:/MS_project/df_lip_example.csv')
# Remove rows where the value in the Molecule column is NaN
df_lip_example = df_lip_example.dropna(subset=['Molecule'])

# 

# Left merge transitions with df on columns "PrecursorName" and "PrecursorAdduct" and "ProductMz"
df_lip_example = df_lip_example.merge(transition_subset, how='left', left_on = ["Molecule", "Precursor Adduct", "Chromatogram Product M/Z"], right_on = ["PrecursorName", "PrecursorAdduct", "ProductMz"])

df_lip_example = df_lip_example.merge(transition_subset, how='left', left_on = ["Precursor Adduct", "Chromatogram Product M/Z"], right_on = ["PrecursorAdduct", "ProductMz"])
df_lip_example = df_lip_example.merge(transitions, how='left', left_on = "Chromatogram Product M/Z", right_on = "ProductMz")


# 