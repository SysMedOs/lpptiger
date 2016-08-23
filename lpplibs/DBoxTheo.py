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

            num_fileds_lst = ['OAP', 'OCP', 'DB', 'OH', 'KETO', 'CHO', 'COOH']

            fa_dct= ox_func(usr_fa_smiles)
            print('fa_dct', fa_dct)
            db_count = fa_dct['DB_count']

            mod_info_df = theodb_oxidizer_cls.mod_sum()
            mod_typ_lst = mod_info_df.columns.tolist()
            mod_sum_df = pd.DataFrame()
            # change the data type from str to int
            mod_info_df.loc['OAP', :] = mod_info_df.loc['OAP', :].astype(int)
            mod_info_df.loc['OCP', :] = mod_info_df.loc['OCP', :].astype(int)
            mod_info_df.loc['DB', :] = mod_info_df.loc['DB', :].astype(int)
            mod_info_df.loc['OH', :] = mod_info_df.loc['OH', :].astype(int)
            mod_info_df.loc['KETO', :] = mod_info_df.loc['KETO', :].astype(int)
            mod_info_df.loc['CHO', :] = mod_info_df.loc['CHO', :].astype(int)
            mod_info_df.loc['COOH', :] = mod_info_df.loc['COOH', :].astype(int)

            # start oxidation
            for db_i in range(1, db_count + 1):
                for _mod in mod_typ_lst:
                    _tmp_mod_lst = mod_sum_df.columns.tolist()
                    if len(_tmp_mod_lst) <= len(mod_typ_lst) and db_i == 1:
                        mod_sum_df[''.join([str(db_i), '-', _mod])] = mod_info_df[_mod]
                    else:
                        for _tmp_mod in _tmp_mod_lst:
                            _tmp_mod_lst = _tmp_mod.split('-')
                            # print('_tmp_mod_lst', _tmp_mod_lst[-2], _tmp_mod_lst)
                            if int(_tmp_mod_lst[-2]) == db_i - 1:
                                _ocp_checker = mod_sum_df.loc['OCP', _tmp_mod]
                                if _ocp_checker == 0:
                                    _new_mod = ''.join([_tmp_mod, '-', str(db_i), '-', _mod])
                                    # print(_new_mod)
                                    # this order is important for the SMILES generated
                                    mod_sum_df.loc[:, _new_mod] = mod_sum_df[_tmp_mod] + mod_info_df[_mod]
                                else:
                                    print('OCP', _ocp_checker, _tmp_mod)
                            else:
                                pass
            mod_sum_df = mod_sum_df.transpose()
            mod_sum_df['MOD_NUM'] = mod_sum_df['OAP'] + mod_sum_df['OCP']

            # filter the OCP and OAP. OAP should be full length
            mod_ocp_sum_df = mod_sum_df.query('OCP == 1')
            # OAP should have all DB, thus MOD_NUM == db_count
            mod_oap_sum_df = mod_sum_df.query('OCP == 0 and MOD_NUM == %d' % db_count)

            # the end of smiles is different for OCP
            mod_ocp_sum_idx_lst = mod_ocp_sum_df.index.tolist()
            mod_ocp_sum_df.loc[mod_ocp_sum_idx_lst, 'FULL_SMILES'] = (fa_dct['DB_pre_part'] +
                                                                      mod_ocp_sum_df['SMILES'] + ')=O')
            mod_ocp_sum_df.is_copy = False
            mod_oap_sum_idx_lst = mod_oap_sum_df.index.tolist()
            mod_oap_sum_df.loc[mod_oap_sum_idx_lst, 'FULL_SMILES'] = (fa_dct['DB_pre_part'] +
                                                                      mod_oap_sum_df['SMILES'] +
                                                                      fa_dct['DB_post_part'])
            mod_oap_sum_df.is_copy = False

            mod_sum_df = mod_ocp_sum_df.append(mod_oap_sum_df)
            _full_smiles_lst = mod_sum_df['FULL_SMILES'].tolist()
            _c_num_lst = []
            for _smiles in _full_smiles_lst:
                _c_num_lst.append(_smiles.count('C'))
            mod_sum_df.loc[:, 'C_NUM'] = _c_num_lst
            mod_sum_df.is_copy = False

            mod_sum_df.to_csv('oxDB.csv')
            num_fileds_lst.append('C_NUM')
            mod_sum_df[num_fileds_lst] = mod_sum_df[num_fileds_lst].astype(str)
            mod_sum_df.loc[:, 'FA_CHECKER'] = (mod_sum_df['C_NUM'] + mod_sum_df['OAP'] + mod_sum_df['OCP'] +
                                               mod_sum_df['DB'] + mod_sum_df['OH'] + mod_sum_df['KETO'] +
                                               mod_sum_df['CHO'] + mod_sum_df['COOH'])
            mod_sum_t_df = mod_sum_df.transpose()
            mod_sum_df.to_csv('oxDB_t.csv')

            print(mod_sum_df.shape)

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


# usr_fa = 'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
usr_fa = 'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'

oxidizer(usr_fa)
