# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

class TheoDB_Oxidizer:
    def __init__(self):
        print("Start to oxidize C=C -->")

    @staticmethod
    def hydroxyl():
        mod_str = r'C(O)C=C'
        pre_mod_counter_dct = {'OH': 1, 'DB': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def hydroxyl_no_db():
        mod_str = r'C(O)CC'
        pre_mod_counter_dct = {'OH': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def keto():
        mod_str = r'C(=O)C=C'
        pre_mod_counter_dct = {'KETO': 1, 'DB': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def keto_no_db():
        mod_str = r'C(=O)CC'
        pre_mod_counter_dct = {'KETO': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def aldehyde():
        mod_str = r'CC=O'
        pre_mod_counter_dct = {'CHO': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def aldehyde_short():
        mod_str = r'C=O'
        pre_mod_counter_dct = {'CHO': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def carboxylic_acid():
        mod_str = r'CC(O)=O'
        pre_mod_counter_dct = {'COOH': 1}
        return mod_str, pre_mod_counter_dct

    @staticmethod
    def carboxylic_acid():
        mod_str = r'C(O)=O'
        pre_mod_counter_dct = {'COOH': 1}
        return mod_str, pre_mod_counter_dct