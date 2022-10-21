import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from sklearn.preprocessing import StandardScaler
import sklearn_json as skljson
import urllib.request
from tqdm import tqdm

from prodec import *


def getCompoundDescriptors(smiles_list : list) :
    """ 
    Compute descriptors for each molecule in the list of smiles.
    Parameters:
        smiles_list: list of smiles
    Returns:
        dictionary of compound descriptors
    """
    descriptors = {}
    print('Computing descriptors for compounds...')
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        desc = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
        descriptors[smiles] = desc

    return descriptors

def getProteinSequences(accession_list : list) :
    """ 
    Get protein sequences for each protein in the list of accessions.
    Parameters:
        accession_list: list of accessions
    Returns:
        dictionary of protein sequences
    """
    sequences = {}
    print('Getting protein sequences...')
    for accession in tqdm(accession_list):
        try:
            url = f'https://www.uniprot.org/uniprot/{accession}.fasta'
            with urllib.request.urlopen(url) as data:
                sequences[accession] = ''.join([line.decode('ascii').strip() for line in data][1:])
        except:
            print('No sequence found for ', accession)

    return sequences

def getProteinDescriptors(accession_list : list, scaler_to_file: str = None, scaler_from_file : str = None):
    """ 
    Compute descriptors for each protein in the list of accessions.
    8 PhysChem descriptors from prodec per amino acid. 
    These descriptors are lineary combined to 50 domains plus the average (401 descriptors in total)
    to make them protein sequence lenght independent.
    A standard scaler is applied to the descriptors.

    Parameters:
        accession_list: list of accessions
        scaler_to_file: path to save the scaler
        scaler_from_file: path to load the scaler
    Returns:
        dictionary of protein descriptors
    """
    
    ndomain = 50
    pdesc = ProteinDescriptors().get_descriptor('PhysChem')
    ndesc = len(pdesc.definition[0])
    pdesc = Transform(TransformType.AVG, pdesc)
    desc = np.zeros((len(accession_list), ndomain*ndesc+1))
    
    sequences = getProteinSequences(accession_list)
    print('Computing descriptors for proteins...')
    for i, accession in enumerate(tqdm(accession_list)):
        desc[i,:] = pdesc.get(sequences[accession], domains=ndomain)

    if scaler_from_file is not None:
        scaler = skljson.from_json(scaler_from_file)
        desc = scaler.transform(desc)
    else:
        scaler = StandardScaler()
        desc = scaler.fit_transform(desc)
    if scaler_to_file is not None:    
        skljson.to_json(scaler, scaler_to_file)

    descriptors = { id : desc[i,:] for i,id in enumerate(accession_list) }

    return descriptors

def getDescriptors( mol_prot_pairs : list, scaler_to_file: str = None, scaler_from_file : str = None):

    """ 
    Compute descriptors for each molecule-target pair in the list of pairs.
    Parameters:
        mol_prot_pairs: list of tuples (smiles, accession)
        scaler_to_file: path to save the scaler
        scaler_from_file: path to load the scaler
    Returns:
        np.array of descriptors of shape (n_pairs, n_descriptors)
    """
    
    # cmp and prot descriptors
    smiles_list = list(set([x[0] for x in mol_prot_pairs]))
    accession_list = list(set([x[1] for x in mol_prot_pairs]))
    cmp_desc = getCompoundDescriptors(smiles_list)
    prot_desc = getProteinDescriptors(accession_list, scaler_to_file, scaler_from_file)
    
    # compound and protein descriptors lenghts
    ncmp = len(cmp_desc[smiles_list[0]])
    nprot = len(prot_desc[accession_list[0]])

    desc = np.zeros((len(mol_prot_pairs), ncmp+nprot))

    for i, (smiles, accession) in enumerate(mol_prot_pairs):
        desc[i,:] = np.concatenate((cmp_desc[smiles], prot_desc[accession]))

    return desc