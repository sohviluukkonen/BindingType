import os
import argparse

import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from .annotation import ClassA_GPCR_HierachicalBindingTypeAnnotation, Kinase_AllostericAnnotation

def add_binding_type_to_papyrus(df : pd.DataFrame, target_type : str = 'GPCR', sources : str = 'both', similarity : bool = False, similarity_th : float = 0.8) -> pd.DataFrame:
    
    """
    Add binding type to Papyrus data frame
    
    Parameters
    ----------
    df : pd.DataFrame
        Papyrus data frame
    target_type : str, optional
        Target type: 'GPCR' or 'Kinase', by default 'GPCR'
    sources : str, optional
        Binding type sources: 'abstracts' (PubMed, PubChem, CrossRef, Patents), 'assays' (ChEMBL) or 'both, by default 'both'
    similarity : bool, optional
        Whether to use similarity to annotate, by default True
    similarity_th : float, optional
        Similarity threshold, by default 0.8

    Returns
    -------
    pd.DataFrame
        Papyrus data frame with binding type
    """

    # Get document and/or assay ids
    if sources not in ['abstracts', 'assays', 'both']:
        raise ValueError('Unknown source type')
    
    if sources == 'both' or sources == 'abstracts':
        document_ids = []
        for doc_ids in df.all_doc_ids.unique():
            for doc_id in doc_ids.split(';'):
                if doc_id.startswith('PMID') or doc_id.startswith('PubChemAID') or doc_id.startswith('DOI') or doc_id.startswith('PATENT'):
                    document_ids.append(doc_id)
        document_ids = list(set(document_ids))
        print('Number of documents:', len(document_ids))

    if sources == 'both' or sources == 'assays':
        chembl_assay_ids = []
        for aids in df.AID.unique():
            for aid in aids.split(';'):
                if aid.startswith('CHEMBL'):
                    chembl_assay_ids.append(aid)
        chembl_assay_ids = list(set(chembl_assay_ids))
        print('Number of ChEMBL assays:', len(chembl_assay_ids))

    # Get binding type by parsing documents and/or assays
    if target_type == 'GPCR':
        parser = ClassA_GPCR_HierachicalBindingTypeAnnotation()
    elif target_type == 'Kinase':
        parser = Kinase_AllostericAnnotation()
    else:
        raise ValueError('Unknown target type')
        
    if sources == 'abstracts':
        document_binding_types = parser(document_ids=document_ids)
        assay_binding_types = None
    elif sources == 'assays':
        assay_binding_types = parser(assay_ids=chembl_assay_ids)
        document_binding_types = None
    elif sources == 'both':
        document_binding_types, assay_binding_types = parser(document_ids=document_ids, assay_ids=chembl_assay_ids)

    # Annotate compounds
    binding_types = list( set(document_binding_types.values()) | set(assay_binding_types.values()) )
    df['BindingType'] = 'Unknown'

    print(f'WARNING: If a compound has multiple binding type annotations, only the first one will be kept.')
    
    for binding_type in binding_types:

        if binding_type == 'Unknown':
            continue

        if document_binding_types is not None:
            for ichikey, dids in zip(df.InChIKey, df.all_doc_ids):    
                
                if any(bd != 'Unknown' for bd in df[df.InChIKey == ichikey].BindingType):
                    continue
                # TODO: implement majority vote
                
                for did in dids.split(';'):
                    if did not in document_binding_types.keys():
                        continue
                    elif document_binding_types[did] == binding_type:
                        df.loc[df.InChIKey == ichikey, 'BindingType'] = binding_type
                        break
        
        if assay_binding_types is not None:
            for ichikey, aids in zip(df.InChIKey, df.AID):    
                
                if any(bd != 'Unknown' for bd in df[df.InChIKey == ichikey].BindingType):
                    continue
                # TODO: implement majority vote / ranking of sources
                
                for aid in aids.split(';'):
                    if aid not in assay_binding_types.keys():
                        continue
                    elif assay_binding_types[aid] == binding_type:
                        df.loc[df.InChIKey == ichikey, 'BindingType'] = binding_type
                        break

    # Annotate compounds with similarity
    if similarity:
        
        # Get canonical SMILES
        df['canonical_SMILES'] = df.SMILES.apply(lambda x: Chem.MolToSmiles(Chem.MolFromSmiles(x), isomericSmiles=False))
        annotated_smiles = df[df.BindingType != 'Unknown'].canonical_SMILES.unique()
        unannotated_smiles = df[df.BindingType == 'Unknown'].canonical_SMILES.unique()

        if len(annotated_smiles) == 0:
            print('WARNING: No annotated compounds found. Similarity annotation will be skipped')
            df.drop('canonical_SMILES', axis=1, inplace=True)
            return df

        # Get fingerprints
        fps_annotated = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 3, nBits=2048) for smi in annotated_smiles]
        fps_unannotated = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 3, nBits=2048) for smi in unannotated_smiles]

        for i, fp in enumerate(fps_unannotated):
            # Get similarities
            similarities = DataStructs.BulkTanimotoSimilarity(fp, fps_annotated)
            
            # If similarity is above threshold, assign binding type to compound with highest similarity to annotated compound  
            if max(similarities) >= similarity_th:
                df.loc[df.canonical_SMILES == unannotated_smiles[i], 'BindingType'] = df[df.canonical_SMILES == annotated_smiles[similarities.index(max(similarities))]].BindingType.values[0]

        df.drop('canonical_SMILES', axis=1, inplace=True)

    return df 
