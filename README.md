
Scripts to annotate GPCR bioactivity binding types and to predict Class A GPCRs binding type based on work in [Burggraaf et al., J. Chem. Inf. Model. (2020)](https://doi.org/10.1021/acs.jcim.0c00695).

# Getting started

## Install

    pip install git+https://github.com/sohviluukkonen/BindingType.git@master

## Binding type annotation

Script to annotate GPCRs bioactivity type based on parsing paper abstracts or PDFs for hierarchical keywords.

Load annotation class: 

    from gpcrBT import *
    annotate = BindingTypeAnnotation()

Parse abstracts from a list of PMIDs:

    doc_ids = [ PMID:15743197, PMID:18307293, PMID:20095623]
    binding_types = annotate.parse_abstracts(document_ids=doc_ids)

Parse abstracts from a .tsv file with PMIDs in `doc_id` columns and save to file .json-file:

    annotate.parse_abstracts(from_file='data/pmids.tsv', out_file='binding_types.json')

Parse PDFs (named PMID:XXXXXX.pdf) in a folder:
    binding_type = annotate.parse_pdfs(input_folder='data/PDFs')


## Class A binding type prediction

Script to predict binary binding type (allo- or orthosteric) of molecule targeting a class A GPCR.

    preds = predictor = BindingTypeClassifier('data/input.tsv')

with `input.tsv` with columns `SMILES` and `accession` containing molecule protein pairs.

# Binding type classifier

A binary proteochemometric Random Forrest-based model was built and trained on Class A GPCR data. 

## Data

Data were obtained by parsing abstracts of molecule-protein pairs of class A GPCRs from CHEMBL22 and filtering out pairs annotated as something else than allo- or orthosteric.

 The model performance was evaluated with a temporal split (train set <= 2013, test set > 2013).

The training set was balanced to correct for the unequal number of orthosteric and allosteric compounds. This was performed using oversampling of the minority class with [SMOTE](https://https://imbalanced-learn.org/stable/references/generated/imblearn.over_sampling.SMOTE.html).

## Descriptors
Molecule features are described with the Morgan bit-fingerprints with a radius of 3 and a length of 1024. The compound descriptors are generated with rdkit.

Proteins features are described by 8 physico-chemical descriptors per residue from [Lenselink et al., J. Cheminform. (2017)](https://doi.org/10.1186/s13321-017-0232-0) that averaged out over 50 domains. The protein descriptors are generated with the [prodec](https://https://github.com/OlivierBeq/ProDEC)-package.


## Model performance on temporal split

| Accuracy | Sensitivity* | Specificity** | MCC  | AUC  |
|----------|--------------|---------------|------|------|
| 0.85     | 0.68         | 0.93          | 0.64 | 0.88 |

*hit rate of allosteric compounds, **hite rate of orthosteric compounds