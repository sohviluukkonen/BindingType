import os
import json
import metapub

import pandas as pd

from io import StringIO
from tqdm.auto import tqdm

from Bio import Entrez
from Bio.Entrez import efetch
Entrez.email = 'A.N.Other@example.com'

from pdfminer.converter import TextConverter
from pdfminer.layout import LAParams
from pdfminer.pdfdocument import PDFDocument
from pdfminer.pdfinterp import PDFResourceManager, PDFPageInterpreter
from pdfminer.pdfpage import PDFPage
from pdfminer.pdfparser import PDFParser
    
def pmid2bibtex(doc_ids, fname):
    
    """ Writes a .bib file from document identification (PMID) in a DataFrame 
    
    Parameters:
        doc_ids (lst) : list of PMID document IDs
        fname (str) : name of ouput .bib file
    """
    
    doc_ids = df.doc_id.unique().tolist()
    fetch = metapub.PubMedFetcher()
    f = open(fname, 'w')
    for doc_id in tqdm(doc_ids):
        try:
            if doc_id.startswith('PMID:'):
                doi = metapub.convert.pmid2doi(doc_id[5:])
                article = fetch.article_by_pmid(doc_id[5:])
                f.write('@article{' + doc_id + ',\n')
                f.write('\tdoi={' + doi + '},\n')
                f.write('\ttitle={' + article.title + '},\n')
                f.write('\tauthor={' + ', '.join(article.authors) + '},\n')
                f.write('\tyear={' + article.year + '},\n')
                f.write('\tjournal={' + article.journal + '},\n')
                f.write('\tabstract={' + doc_id + '},\n')
                f.write('}\n')
        except:
            # If doc_id is not a PMID, pass 
            print(doc_id)
    f.close()

class BindingTypeAnnotation():
    
    """ Annotate binding type of a GPCR publications by parsing abstracts or PDFs for hierarchy of keywords """

    def abstract_parser(self, document_ids : list = None, input_file : str = None, output_file : str = None):

        """ Parse abstracts for binding type keywords
        
        Parameters:
            document_ids (list of str) : list of document IDs (PMID)
            input_file (str) : path to input file
            output_file (str) : path to output file
        """
        
        if input_file is not None:
            document_ids = pd.read_csv(input_file, sep='\t').doc_id.unique().tolist()
        elif document_ids is None:
            raise ValueError('No document IDs provided. Provide either a list of document IDs or a path to a file containing document IDs.')

        # Get, parse and annotate abstracts
        dct_binding_types = {}
        for pmid in tqdm(document_ids, desc='Parsing abstracts'):
            try :
                handle = efetch(db='pubmed', id=pmid[5:], retmode='text', rettype='abstract')
                abstract = handle.read().lower()
                dct_binding_types[pmid] = self.hierarchichal_binding_site_attribution(abstract)
            except:
                print('Could parse {}'.format(pmid))

        # Write to file
        if output_file is not None:
            with open(output_file, 'w') as f: json.dump(dct_binding_types, f)
        
        return dct_binding_types

    def pdf_parser(self, input_folder : str = None, output_file : str = None):

        """ Parse PDFs for binding type keywords
        
        Parameters:
            input_folder (str) : path to input folder
            output_file (str) : path to output file
        """
        
        if input_folder is None:
            raise ValueError('No input folder provided. Provide a path to a folder containing PDFs.')
        files = [f for f in os.listdir(input_folder) if f.endswith('.pdf')]

        # Get, parse and annotate PDFs
        dct_binding_types = {}
        for pmid in tqdm(os.listdir(files), desc='Parsing PDFs'):
            text = self.convert_pdf_to_string(os.path.join(input_folder, fname))
            text = text.lower().split("references",1)[0]
            dct_binding_types[pmid] = self.hierarchichal_binding_site_attribution(text)
        
        # Write to file
        if output_file is not None:
            with open(output_file, 'w') as f: json.dump(dct_binding_types, f)
        
        return dct_binding_types

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
        keywords = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'keywords.json'), 'r'))
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

    def convert_pdf_to_string(self, file_path):
        
        """ Converts a PDF file to string
        
        Parameters:
            file_path (str): path to PDF file
        Returns:
            (str) : str containing text from PDF"""

        output_string = StringIO()
        with open(file_path, 'rb') as in_file:
            parser = PDFParser(in_file)
            doc = PDFDocument(parser)
            rsrcmgr = PDFResourceManager()
            device = TextConverter(rsrcmgr, output_string, laparams=LAParams())
            interpreter = PDFPageInterpreter(rsrcmgr, device)
            for page in PDFPage.create_pages(doc):
                interpreter.process_page(page)

        return(output_string.getvalue())
    