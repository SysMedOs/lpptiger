# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re

import numpy as np
import pandas as pd

# from rdkit import Chem
# from rdkit.Chem import Draw


class TheoDB_Oxidizer:
    def __init__(self):
        print("Start to oxidize C=C -->")

    @staticmethod
    def mod_sum():
        mod_info_df = pd.read_csv('ModConfig.csv', index_col=0, dtype={'OAP': np.int32})
        # print(mod_info_df)
        return mod_info_df

# construct a decorator
def bulk_oxidizer(theodb_oxidizer_cls):
    """

    :param theodb_cls:
    :return:
    """

    def _bulk_oxidizer(ox_func):
        def __bulk_oxidizer(usr_fa_smiles):

            mod_str = ''

            fa_dct= ox_func(usr_fa_smiles)
            print('fa_dct', fa_dct)
            db_count = fa_dct['DB_count']

            mod_info_df = theodb_oxidizer_cls.mod_sum()
            mod_typ_lst = mod_info_df.columns.tolist()
            mod_counter_dct = {}
            mod_sum_df = pd.DataFrame()
            mod_info_df.loc['OAP', :] = mod_info_df.loc['OAP', :].astype(int)
            mod_info_df.loc['OCP', :] = mod_info_df.loc['OCP', :].astype(int)
            mod_info_df.loc['DB', :] = mod_info_df.loc['DB', :].astype(int)
            mod_info_df.loc['OH', :] = mod_info_df.loc['OH', :].astype(int)
            mod_info_df.loc['KETO', :] = mod_info_df.loc['KETO', :].astype(int)
            mod_info_df.loc['CHO', :] = mod_info_df.loc['CHO', :].astype(int)
            mod_info_df.loc['COOH', :] = mod_info_df.loc['COOH', :].astype(int)

            for db_i in range(1, db_count + 1):
                for _mod in mod_typ_lst:
                    _tmp_mod_lst = mod_sum_df.columns.tolist()
                    if len(_tmp_mod_lst) <= len(mod_typ_lst) and db_i == 1:
                        mod_sum_df[''.join([str(db_i), '_', _mod])] = mod_info_df[_mod]
                    else:
                        for _tmp_mod in _tmp_mod_lst:
                            if int(_tmp_mod[0]) < db_i:
                                _new_mod = ''.join([str(db_i), '_', _mod, '_', _tmp_mod])
                                print(_new_mod)
                                mod_sum_df[_new_mod] = mod_info_df[_mod] + mod_sum_df[_tmp_mod]

            print(mod_typ_lst)
            print(mod_sum_df)


        return __bulk_oxidizer
    return _bulk_oxidizer


@bulk_oxidizer(TheoDB_Oxidizer)
def oxidizer(usr_fa_smiles):

    db_info_dct = {}

    fa_rgx = re.compile(r'(OC[(])([C/\\=()]*)(C[)]=O)')

    db_rgx = re.compile(r'(C/C[=]C\\)')

    # Construct regular expression for DB
    # if use FA: 'OC(CCCCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
    # get ('CCCCCC', 'C/C=C\\C/C=C\\C/C=C\\C/C=C\\', 'C/C=C\\', 'CCCC')
    fa_main_rgx = re.compile(r'(C*)((C[/]C[=]C\\){1,8})(C*)')

    fa_checker = re.match(fa_rgx, usr_fa_smiles)
    if fa_checker:
        pre_fa_lst = fa_checker.groups()
        fa_main_str = pre_fa_lst[1]
        fa_main_checker = re.match(fa_main_rgx, fa_main_str)
        if fa_main_checker:
            pre_db_lst = fa_main_checker.groups()

            if len(pre_db_lst) >= 3:
                db_str = pre_db_lst[1]
                db_counter = db_str.count('C/C=C\\')
                print('db_str', db_str, 'db_counter', db_counter)

                db_info_dct['DB_count'] = db_counter
                db_info_dct['DB_pre_part'] = ''.join(['OC(', pre_db_lst[0]])
                db_info_dct['DB_main_part'] = pre_db_lst[1]
                db_info_dct['DB_post_part'] = ''.join([pre_db_lst[-1], 'C)=O'])
                if (''.join([db_info_dct['DB_pre_part'], db_info_dct['DB_main_part'], db_info_dct['DB_post_part']])
                        == usr_fa_smiles):
                    db_info_dct['DB_full_fa'] = usr_fa_smiles
                    db_info_dct['DB_C_count'] = usr_fa_smiles.count('C')
                    db_info_dct['DB_start'] = db_info_dct['DB_pre_part'].count('C') + 2  # for that 2 C in DB pattern
                    db_info_dct['DB_omega'] = db_info_dct['DB_post_part'].count('C') + 1  # for that C in DB pattern

            else:
                pass
        else:
            pass

        return db_info_dct


usr_fa = 'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'

oxidizer(usr_fa)
