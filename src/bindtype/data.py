import pandas as pd
from tqdm import tqdm

# Mapping of ChEMBL ID to accessions and protein classification
data = pd.read_table('../data/preprocessingFiles/S5_Compound_Dataset.tsv.gz')
chembl22_id2accession = pd.read_table('preprocessingFiles/Chembl22_id2accession.tsv')
chembl22_protein_class = pd.read_table('preprocessingFiles/Chembl22_protein_classification.tsv')

for chembl_id in tqdm(data.CHEMBL_ID_protein.unique()):
    data.loc[data.CHEMBL_ID_protein == chembl_id, 'accession'] = chembl22_id2accession.loc[chembl22_id2accession.chembl_id == chembl_id, 'accession'].values[0]
    for i in range(1,9):
        data.loc[data.CHEMBL_ID_protein == chembl_id, f'L{i}'] = chembl22_protein_class.loc[chembl22_protein_class.chembl_id == chembl_id, f'l{i}'].values[0]

# Fix obsolete accession code
data.loc[data.accession == 'P10845', 'accession'] = 'P0DPI0'

# Save data
data = data[['CHEMBL_ID_compound', 'canonical_smiles', 'CHEMBL_ID_protein', 'accession', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'year', 'binding_type', 'activity_id', 'pchembl_value']]
data.to_csv('../data/BindingTypeDataset.tsv', sep='\t', index=False)

# Keep only class A GPCR
print('Number of data points:', len(data))
data = data[data.L2 == 'Family A G protein-coupled receptor']
print('Number of data points after filtering for class A GPCRs:', len(data))

# Keep only data points with known binding type
print('Number of data points: ', len(data))
data = data[data.binding_type.isin(['allosteric', 'orthosteric'])].reset_index(drop=True)
print('Number of data points after filtering for allo- and orthosteric compounds: ', len(data))

# Binary binding type
data['binding_binary'] = data['binding_type'].apply(lambda x: 1 if x == 'allosteric' else 0)

# Temporal split without data-leakage
lim_year = 2013
data.loc[ data.year <= lim_year, 'Subset'] = 'train'
data.loc[ data.year > lim_year, 'Subset'] = 'test'

duplicates = data[data.CHEMBL_ID_compound.duplicated()].CHEMBL_ID_compound.unique()
for cmp in duplicates:
    years = data[data.CHEMBL_ID_compound == cmp].year.tolist()
    if min(years) <= lim_year : data.loc[data.CHEMBL_ID_compound == cmp, 'Subset'] = 'train'
    else: data.loc[data.CHEMBL_ID_compound == cmp, 'Subset'] = 'test'

# Save data
data.to_csv('../data/BindingTypeDataset_ClassAGPCRs.tsv', sep='\t', index=False)