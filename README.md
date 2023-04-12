
Scripts to annotate GPCR bioactivity binding types and to predict Class A GPCRs binding type based on work in [Burggraaf et al., J. Chem. Inf. Model. (2020)](https://doi.org/10.1021/acs.jcim.0c00695).

# Getting started

## Install

    pip install git+https://github.com/sohviluukkonen/BindingType.git@main

## Binding type annotation

Script to annotate GPCRs bioactivity type based on parsing paper abstracts or PDFs for hierarchical keywords.

Load annotation class: 

    from bindtype import *
    annotate = BindingTypeAnnotation()

Parse abstracts from a list of PMIDs:

    doc_ids = [ PMID:15743197, PMID:18307293, PMID:20095623]
    binding_types = annotate.parse_abstracts(document_ids=doc_ids)

Parse abstracts from a .tsv file with PMIDs in `doc_id` columns and save to file .json-file:

    annotate.parse_abstracts(from_file='pmids.tsv', out_file='binding_types.json')

Parse PDFs (named PMID:XXXXXX.pdf) in a folder:
    binding_type = annotate.parse_pdfs(input_folder='PDF_folder')


## Class A binding type prediction

The PCM model for binding type predictions can be found in [Classifier](/github/BindingType/tree/Classifier) branch.