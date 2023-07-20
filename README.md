
Package to annote binding type of bioactivity measures based on keyword search of
- __abstracts__ from PubMed, PubChem assay description, CrossRef or Google Patents
- __assay descriptions__ from ChEMBL assay descriptions

The annotation is currently supported for two types of targets
- __Class A GPCRs__ with a 3-level hierarchical keyword search to annotate compounds as _orthosteric, allosteric, bitopic, covalent_ or  _unknown_ based work in [Burggraaf et al., J. Chem. Inf. Model. (2020)](https://doi.org/10.1021/acs.jcim.0c00695)
- __Protein Kinases__ with 1-level keyword search to annotate compounds as _allosteric_ or _unknown_ based (extended) keywords from [Christmann-Franck et al., J. Chem. Inf. Model. (2016)](https://doi.org/10.1021/acs.jcim.6b00122)

# Getting started

## Install

    pip install git+https://github.com/sohviluukkonen/BindingType.git@main

## Usage

The package has both an API and a CLI which can process either
- [Papyrus](https://github.com/OlivierBeq/Papyrus-scripts) datasets
- lists of document and/or assay IDs

### Papyrus data

In the case of Papyrus-dataframe, the annotation will a new `BindingType` column to the dataframe and can be done from the command line with
```
bindtype_papyrus -i <dataset.csv/.tsv> -tt <GPCR/Kinase>
```
or with the API with
```
from bindtype.papyrus import add_binding_type_to_papyrus
df = add_binding_type_to_papyrus(df, target_type=GPCR/Kinase)
```

There is also an option to annotate all 'unknown' compounds that based on their Tanimoto similarity to the annotated compounds: `-sim, --similarity` flag in the CLI and `similarity=True` in the API.

### General usage

In the more general case, the annotation will create dictionaries based list of document IDs and/or assays IDs. This can be done either from the command line with
```
bindtype -did <document_id_file_path> -aid <assay_id_file_path> -tt <GPCR/Kinase>
```
or with the API with 
```
# for the GPCRs
from bindtype import ClassA_GPCR_HierachicalBindingTypeAnnotation
parser = ClassA_GPCR_HierachicalBindingTypeAnnotation()

# for the kinases
from bindtype import Kinase_AllostericAnnotation
parser = Kinase_AllostericAnnotation()

# Only abstracts
dct_doc_annotations = parser(document_ids=list_of_document_ids)

# Only assay descriptions
dct_assay_annotations = parser(assay_ids=list_of_assay_ids)

# Both
dct_doc_annotations, dct_assay_annotations = parser(document_ids=list_of_document_ids, assay_ids=list_of_assay_ids)
```

As the scripts were developed with data from [Papyrus](https://github.com/OlivierBeq/Papyrus-scripts) and uses document and assay description IDs should be in the format used in the `all_doc_ids` and `AID` columns: __PMID:<pubchem_id>__, __PubChemAID:<pubchem_assay_id>__, __DOI:<doi>__, __PATENT:<patent_id>__ and __<chembl_assay_id>__.
