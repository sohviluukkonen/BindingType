import os
import json
import argparse
import requests

from bs4 import BeautifulSoup
from tqdm.auto import tqdm
from typing import List, Dict

from Bio import Entrez
from Bio.Entrez import efetch
Entrez.email = 'A.N.Other@example.com'

from chembl_webresource_client.new_client import new_client


class BindingTypeAnnotation():
    
    """ Annotate binding type of a based on keywords in abstracts or assay descriptions """

    def __init__(self):
        self.documents = {}
        self.assays = {}

    def retrieve_abstracts_from_PubMed(self):

        """ Retrieve abstracts from PubMed using document IDs (PMID) """

        pmids = [doc_id for doc_id in self.document_ids if doc_id.startswith('PMID:')]
                 
        if len(pmids) == 0:
            print('No PubMed IDs provided. Parsing PubMed will be skipped.')
        else:
            for pmid in tqdm(pmids, desc='Parsing PubMed'):
                try:
                    handle = efetch(db='pubmed', id=self.pmid[5:], retmode='text', rettype='abstract')
                    self.documents[pmid] = handle.read().lower()
                except:
                    print('Could not retrieve abstract for {}'.format(pmid))

    def retrieve_descriptions_from_PubChem(self):

        """ Retrieve assay descriptions from PubChem using document IDs (PubChemAID) """

        aids = [doc_id for doc_id in self.document_ids if doc_id.startswith('PubChemAID:')]
                 
        if len(aids) == 0:
            print('No PubChemAID IDs provided. Parsing PubChem will be skipped.')
        else:
            for aid in tqdm(aids, desc='Parsing PubChem'):
                try:
                    handle = efetch(db='pcassay', id=self.aid[11:], retmode='text', rettype='abstract')
                    self.documents[aid] = handle.read().lower()
                except:
                    print('Could not retrieve description for {}'.format(aid))

    def retrieve_abstracts_from_CrossRef(self):

        """ Retrieve abstracts from CrossRef using document IDs (DOI) """

        dois = [doc_id for doc_id in self.document_ids if doc_id.startswith('DOI:')]
                 
        if len(dois) == 0:
            print('No DOI IDs provided. Parsing CrossRef will be skipped.')
        else:
            for doi in tqdm(dois, desc='Parsing CrossRef'):
                try:
                    url = f'https://api.crossref.org/works/{doi[4:]}'
                    r = requests.get(url)
                    crossref = r.json()
                    self.documents[doi] = crossref['message']['abstract'].lower()
                except:
                    print('Could not retrieve abstract for {}'.format(doi))
                
    def retrieve_abstracts_from_GooglePatents(self):

        """ Retrieve abstracts from Google Patents using document IDs (Patent number) """

        patents = [doc_id for doc_id in self.document_ids if doc_id.startswith('PATENT:')]
                 
        if len(patents) == 0:
            print('No patent IDs provided. Parsing Google Patents will be skipped.')
        else:
            for patent in tqdm(patents, desc='Parsing Google Patents'):
                try:
                    url = f'https://patents.google.com/patent/{patent[7:]}/en'
                    r = requests.get(url)
                    soup = BeautifulSoup(r.text, 'html.parser')
                    meta = soup.find_all("meta")

                    for m in meta:
                        if m.get("name") == "description":
                           self.documents[patent] = m.get("content").lower()
                           break
                except:
                    print('Could not retrieve abstract for {}'.format(patent))


    def retrieve_description_from_ChEMBL(self):

        """ Retrieve assay description from ChEMBL using document IDs (ChEMBL ID) """

        chembl_aids = [assay_id for assay_id in self.assay_ids if assay_id.startswith('CHEMBL:')]
                 
        if len(chembl_aids) == 0:
            print('No ChEMBL IDs provided. Parsing ChEMBL will be skipped.')
        else:
            chembl_assays = new_client.assay
            descriptions = chembl_assays.filter(assay_id__in=chembl_aids).only(['description'])
            for chembl_aid, description in tqdm(zip(chembl_aids, descriptions), total=len(chembl_aids), desc='Parsing ChEMBL assay descriptions'):
                description = description['description']
                if description is not None:
                    self.assays[chembl_aid] = description.lower()
                else:
                    print('Could not retrieve description for {}'.format(chembl_aid))

    def run_in_parallel(self, *funcs):
        """ Run functions in parallel """

        from multiprocessing import Pool
        from functools import partial

        with Pool() as p:
            p.map(partial(funcs))

    def drop_unknowns(self, d : Dict) -> Dict:
        """ Drop unknown binding types """

        return {k: v for k, v in d.items() if v != 'Unknown'}
    
    def __call__(self, document_ids: List[str] = None, assay_ids: List[str] = None, out_prefix: str = None, keep_unknowns: bool = False):

        if document_ids is None and assay_ids is None:
            raise ValueError('Provide either document IDs or ChEMBL assay IDs.')

        if document_ids:
            self.document_ids = document_ids

            self.retrieve_abstracts_from_PubMed()
            self.retrieve_abstracts_from_CrossRef()
            self.retrieve_abstracts_from_GooglePatents()
            self.retrieve_descriptions_from_PubChem()

            document_binding_types = {}
            for doc_id, text in self.documents.items():
                document_binding_types[doc_id] = self.annotate(text)

            if not keep_unknowns:
                document_binding_types = self.drop_unknowns(document_binding_types)

            if out_prefix is not None:
                with open(out_prefix + '_docs.json', 'w') as f:
                    json.dump(document_binding_types, f)
        
        if assay_ids:
            self.assay_ids = assay_ids

            self.retrieve_description_from_ChEMBL()

            assay_binding_types = {}
            for assay_id, text in self.assays.items():
                assay_binding_types[assay_id] = self.annotate(text)

            if not keep_unknowns:
                assay_binding_types = self.drop_unknowns(assay_binding_types)

            if out_prefix is not None:
                with open(out_prefix + '_assays.json', 'w') as f:
                    json.dump(assay_binding_types, f)

        if document_ids and assay_ids:
            return document_binding_types, assay_binding_types
        elif document_ids:
            return document_binding_types
        elif assay_ids:
            return assay_binding_types


class Kinase_AllostericAnnotation(BindingTypeAnnotation):

    """ Annotate binding type of kinases based on keywords in abstracts or assay descriptions """

    def __init__(self):
        super(Kinase_AllostericAnnotation).__init__()
        self.annotate = self.annotate_allosteric

    def annotate_allosteric(self, text: str) -> str:

        """ Annotate binding type of kinases based on keywords in abstracts or assay descriptions

        Parameters
        ----------
        text : str
            Abstract or assay description

        Returns
        -------
        str
            Binding type
        """

        text = text.lower()

        keywords = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Kinase_allosteric_keywords.json'), 'r'))
        for k, v in keywords.items(): keywords[k] = [item.lower() for item in v]

        if any([keyword in text for keyword in keywords]):
            return 'Allosteric'
        else:
            return 'Unknown'


class ClassA_GPCR_HierachicalBindingTypeAnnotation(BindingTypeAnnotation):
    
    """ Annotate binding type of Class A GPCR ligands based on hierarchical keywords in abstracts or assay descriptions """

    def __init__(self):
        super(ClassA_GPCR_HierachicalBindingTypeAnnotation).__init__()
        self.annotate = self.hierarchichal_binding_site_attribution


    def hierarchichal_binding_site_attribution(self, text):
        
        """ Attributes compound binding type by text mining a string with hierarchical classification on keywords
        
        Parameters:
            text (str): str to be mined
        
        Returns:
            site (str): 'Orthosteric', 'Allosteric', 'Bitopic', 'Covalent' or 'Unknown'
        """
        text = text.lower()
        site = None

        # Load keywords
        keywords = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ClassA_GPCR_hierarchical_keywords.json'), 'r'))
        for k, v in keywords.items(): keywords[k] = [item.lower() for item in v]
        
        if any(item in text for item in keywords['Bitopic']): # Check if ligand binds to both sites
            site = 'Bitopic'
        
        elif any(item in text for item in keywords['Covalent']): # Check if ligand binds covalently
            site = 'Covalent'
        
        else: # If not, -- 1. level allo/ortho parsing
            
            if any(item in text for item in keywords['Allo1']): # Contains 1. level allosteric keyword
                
                if not any(item in text for item in keywords['Ortho1']): # Does not contain a 1. level orthosteric keyword
                    site = 'Allosteric'
                else: # Contains both types of 1. level keywords
                    x = self.compare_frequency_in_string(text, keywords['Allo1'], keywords['Ortho1'])
                    if x == 0 : site = 'Allosteric'
                    if x == 1 : site = 'Orthosteric'

            elif any(item in text for item in keywords['Ortho1']): # Contains only 1. level orthosteric keyword
                site = 'Orthosteric'
            
        if site == None: # Does not contain any or same number 1. level keywords -- 2. level allo/ortho parsing

            if any(item in text for item in keywords['Allo2']): # Contains 2. level allosteric keyword

                if not any(item in text for item in keywords['Ortho2']): # Does not contain a 2. level orthosteric keyword
                    site = 'Allosteric'
                else: # Contains both types of 2. level keywords
                    x = self.compare_frequency_in_string(text, keywords['Allo2'], keywords['Ortho2'])
                    if x == 0 : site = 'Allosteric'
                    if x == 1 : site = 'Orthosteric'

            elif any(item in text for item in keywords['Ortho2']): # Contains only 2. level orthosteric keyword
                site = 'Orthosteric'

        if site == None: # Does not contain any or same number 2. level keywords -- 3. level allo/ortho parsing
            
            if any(item in text for item in keywords['Ortho3']):
                site = 'Orthosteric'
            else:
                site = 'Unknown'
                
        return site

    def compare_frequency_in_string(self, text, lst1, lst2):
        
        """ Counts and compares how times elements of two lists of strings are presenet in a string.
        
        Paremeters:
            text (str)
            lst1 (list of str)
            lst2 (list of str)
        Returns:
            0 if elements of lst1 more present than of lst2 in text
            1 if elements of lst2 more present that of lst1 in text
            2 if elements of lst1 and lst2 equally present in text
        """
        
        count1, count2 = 0, 0
        for item in lst1: count1 += text.count(item)
        for item in lst2: count2 += text.count(item)
        if count1 > count2: return 0
        elif count1 < count2: return 1
        else: return 2    



