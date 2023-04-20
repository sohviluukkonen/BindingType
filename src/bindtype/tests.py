import os
import pandas as pd

from unittest import TestCase

from .annotation import ClassA_GPCR_HierachicalBindingTypeAnnotation, Kinase_AllostericAnnotation
from .papyrus import add_binding_type_to_papyrus

# TODO : add test for differnt methods in annotation.py
# class TestAnnotation(TestCase):

class TestPapyrus(TestCase):

    def setUp(self):
        self.df_AR = pd.read_csv(os.path.join(os.path.dirname(__file__), 'test_data_AR.csv'))

    def test_kinase_with_similarity(self):
        df = add_binding_type_to_papyrus(self.df_AR, target_type='Kinase', sources='both', similarity=True, similarity_th=0.8)

    def test_kinase_without_similarity(self):
        df = add_binding_type_to_papyrus(self.df_AR, target_type='Kinase', sources='both', similarity=False)

    def test_gpcr_with_similarity(self):
        df = add_binding_type_to_papyrus(self.df_AR, target_type='GPCR', sources='both', similarity=True, similarity_th=0.8)

    def test_gpcr_without_similarity(self):
        df = add_binding_type_to_papyrus(self.df_AR, target_type='GPCR', sources='both', similarity=False)
        



        

    