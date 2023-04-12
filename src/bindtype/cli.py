import argparse

from .annotation import ClassA_GPCR_HierachicalBindingTypeAnnotation, Kinase_AllostericAnnotation
from .papyrus import add_binding_type_to_papyrus

def bindtype_cli_json():

    """
    Command line interface for binding type annotation of documents. 
    The annotated document IDs are saved in a JSON file.
    """
    
    parser = argparse.ArgumentParser(description='Annotate binding type of kinases based on keywords in abstracts or assay descriptions')
    parser.add_argument('-dip', '--document_ids_path', type=str, default=None, help='Path to file with document IDs (PMID, PubChemAID, DOI, PATENT)')
    parser.add_argument('-aip', '--assay_ids_path', type=str, default=None, help='Path to file with assay IDs (ChEMBL ID)')
    parser.add_argument('-tt', '--target_type', type=str, default='GPCR', help="Target type ('Kinase' or 'GPCR')")
    parser.add_argument('-o', '--out_prefix', type=str, help='Output file prefix')
    parser.add_argument('--keep_unknowns', action='store_true', help='Keep unknown binding types')
    args = parser.parse_args()

    if args.document_ids_path is None and args.assay_ids_path is None:
        raise ValueError('File with either document IDs or assay IDs must be provided')
    
    if args.target_type == 'GPCR':
        parser = ClassA_GPCR_HierachicalBindingTypeAnnotation()
    elif args.target_type == 'Kinase':
        parser = Kinase_AllostericAnnotation()
    else:
        raise ValueError('Target type must be either "GPCR" or "Kinase"')
    
    if args.document_ids_path and args.assay_ids_path:
        _, _ = parser(document_ids_path=args.document_ids_path, assay_ids_path=args.assay_ids_path, out_prefix=args.out_prefix, keep_unknowns=args.keep_unknowns)
    elif args.document_ids_path:
        _ = parser(document_ids_path=args.document_ids_path, out_prefix=args.out_prefix, keep_unknowns=args.keep_unknowns)
    elif args.assay_ids_path:
        _ = parser(assay_ids_path=args.assay_ids_path, out_prefix=args.out_prefix, keep_unknowns=args.keep_unknowns)

def bindtype_cli_papyrus():

    """
    Command line interface for binding type annotation of Papyrus data.
    The input file should be a csv or tsv file in the 'Papyrus' format with the following columns:
    - SMILES
    - InChIKey
    - all_doc_ids
    - all_assay_ids
    The output file is a csv or tsv file with the same columns as the input file and an additional column 'BindingType'.
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True, help='File containing Papyrus input data')
    parser.add_argument('-o', '--output', type=str, required=False, help='File to save Papyrus data with binding type')
    parser.add_argument('-t', '--target_type', type=str, required=False, default='GPCR', help='Target type (GPCR or Kinase)')
    parser.add_argument('-s', '--sources', type=str, required=False, default='both', help='Sources to use for annotation (abstracts, assays, both)')
    parser.add_argument('-sim', '--similarity', type=bool, required=False, default=False, help='Use similarity to annotate compounds')
    parser.add_argument('-sim_th', '--similarity_threshold', type=float, required=False, default=0.8, help='Similarity threshold')
    args = parser.parse_args()

    # Load data
    if args.input.endswith('.csv'):
        df = pd.read_csv(args.input)
        suffix = '.csv'
    elif args.input.endswith('.tsv'):
        df = pd.read_table(args.input)
        suffix = '.tsv'
    else:
        raise ValueError('Unknown file format: should be .csv or .tsv')

    # Add binding type
    df = add_binding_type_to_papyrus(df, target_type=args.target_type, sources=args.source, similarity=args.similarity, similarity_th=args.similarity_threshold)

    # Save data
    if args.output is None:
        df.to_csv(args.input.replace(suffix, '_with_binding_type' + suffix), index=False, sep='\t' if suffix == '.tsv' else ',')
    else:
        df.to_csv(args.output, index=False, sep='\t' if args.output.endswith('.tsv') else ',') 
