# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re
import json

import numpy as np
import pandas as pd

from AbbrGenerator import AbbrGenerator

# from rdkit import Chem
# from rdkit.Chem import Draw


class TheoDB_Oxidizer:
    def __init__(self):
        print("Start to oxidize C=C -->")

    @staticmethod
    def mod_sum(usr_mod_table):
        mod_info_df = pd.read_csv(usr_mod_table, index_col=0, dtype={'OAP': np.int32})
        # print(mod_info_df)
        return mod_info_df

# construct a decorator
def bulk_oxidizer(theodb_oxidizer_cls):
    """

    :param theodb_cls:
    :return:
    """

    def _bulk_oxidizer(ox_func):
        def __bulk_oxidizer(usr_fa_dct, usr_mod_table):

            nam_fields_last = ['OAP', 'OCP', 'DB', 'OH', 'KETO', 'CHO', 'COOH']

            fa_dct = ox_func(usr_fa_dct)
            # print('fa_dct', fa_dct)
            db_count = fa_dct['DB_count']
            if fa_dct['DB_LINK_type'] in ['O-', 'P-']:
                ocp_end_part = ''
            else:
                ocp_end_part = ')=O'

            mod_info_df = theodb_oxidizer_cls.mod_sum(usr_mod_table)
            mod_typ_lst = mod_info_df.columns.tolist()
            mod_sum_df = pd.DataFrame()
            
            if db_count > 0:
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
                                        # print('OCP', _ocp_checker, _tmp_mod)
                                        pass
                                else:
                                    pass
                mod_sum_df = mod_sum_df.transpose()
                # no more than one keto & max reduce 1 DB
                _min_DB = db_count - 1
                mod_sum_df = mod_sum_df[(mod_sum_df.KETO < 2) & (mod_sum_df.DB >= _min_DB)]
                mod_sum_df['MOD_NUM'] = mod_sum_df['OAP'] + mod_sum_df['OCP']
    
                # filter the OCP and OAP. OAP should be full length
                mod_ocp_sum_df = mod_sum_df.query('OCP == 1')
                # OAP should have all DB, thus MOD_NUM == db_count
                mod_oap_sum_df = mod_sum_df.query('OCP == 0 and MOD_NUM == %d' % db_count)
    
                # the end of smiles is different for OCP
                mod_ocp_sum_idx_lst = mod_ocp_sum_df.index.tolist()
                mod_ocp_sum_df.loc[mod_ocp_sum_idx_lst, 'FULL_SMILES'] = (fa_dct['DB_pre_part'] +
                                                                          mod_ocp_sum_df['SMILES'] + ocp_end_part)
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
    
                # mod_sum_df.to_csv('oxDB.csv')
                nam_fields_last.append('C_NUM')
                mod_sum_df[nam_fields_last] = mod_sum_df[nam_fields_last].astype(str)

                # here CHO@C & COOH@C is 1 or 0 for T/F
                mod_sum_df.loc[:, 'FA_CHECKER'] = (fa_dct['DB_LINK_type'] + mod_sum_df['C_NUM'] + ':'
                                                   + mod_sum_df['DB'] + '[' + mod_sum_df['DB'] + 'xDB,'
                                                   + mod_sum_df['OH'] + 'xOH,'
                                                   + mod_sum_df['KETO'] + 'xKETO]<CHO@C'
                                                   + mod_sum_df['CHO'] + ',COOH@C'
                                                   + mod_sum_df['COOH'] + '>{OAP:'
                                                   + mod_sum_df['OAP'] + ',OCP:'
                                                   + mod_sum_df['OCP'] + '}')
                abbr_gen = AbbrGenerator()
                mod_sum_df['FA_ABBR'] = ''
                mod_sum_df['FA_TYPE'] = ''
                mod_sum_df['FA_JSON'] = ''
                for (_fa_idx, _fa_row) in mod_sum_df.iterrows():
                    _fa_code = str(_fa_row['FA_CHECKER'])
                    _fa_abbr, _fa_typ = abbr_gen.decode(_fa_code)
                    _fa_row['FA_ABBR'] = _fa_abbr
                    # print('_fa_code', _fa_code)
                    # print(_fa_typ, ' | ', _fa_abbr)
                    _fa_row['FA_TYPE'] = _fa_typ

                    _fa_row['FA_JSON'] = json.dumps({'LINK_TYPE': fa_dct['DB_LINK_type'], 'C': _fa_row['C_NUM'],
                                                     'DB': _fa_row['DB'], 'OH': _fa_row['OH'],
                                                     'KETO': _fa_row['KETO'], 'CHO': _fa_row['CHO'],
                                                     'COOH': _fa_row['COOH'], 'OAP': _fa_row['OAP'],
                                                     'OCP': _fa_row['OCP']})
                    # print(_fa_row['FA_JSON'])
                # mod_sum_t_df = mod_sum_df.transpose()
                # print(mod_sum_t_df.columns.tolist())
                # mod_sum_df.to_csv('oxDB_t.csv')

            else:
                if fa_dct['DB_LINK_type'] == 'O-':
                    _lyso_fa_abbr = 'O-%i:0' % fa_dct['DB_C_count']
                elif fa_dct['DB_LINK_type'] == 'P-':
                    _lyso_fa_abbr = 'P-%i:0' % fa_dct['DB_C_count']
                elif fa_dct['DB_LINK_type'] == '':
                    _lyso_fa_abbr = '%i:0' % fa_dct['DB_C_count']
                else:
                    _lyso_fa_abbr = '%i:0' % fa_dct['DB_C_count']

                unox_json = ('{"C": "%i", "KETO": "0", "OH": "0", "OAP": "0", "OCP": "0", "COOH": "0", "DB": "0", '
                             '"LINK_TYPE": "", "CHO": "0"}' % fa_dct['DB_C_count'])

                unox_dct = {'SMILES': fa_dct['DB_full_fa'], 'OAP': 0, 'OCP': 0, 'DB': 0,
                            'OH': 0, 'KETO': 0, 'CHO': 0, 'COOH': 0, 'MOD_NUM': 0,
                            'FULL_SMILES': fa_dct['DB_full_fa'], 'C_NUM': fa_dct['DB_C_count'],
                            'FA_CHECKER': '%i:0[0xDB,0xOH,0xKETO]<CHO@C0,COOH@C0>{OAP:0,OCP:0}' % fa_dct['DB_C_count'],
                            'FA_ABBR': _lyso_fa_abbr, 'FA_TYPE': 'UNMOD', 'FA_JSON': unox_json}

                mod_sum_df = pd.DataFrame(unox_dct, index=['0-no_oxidation'])

            lyso_json = ('{"C": "0", "KETO": "0", "OH": "0", "OAP": "0", "OCP": "0", "COOH": "0", "DB": "0", '
                         '"LINK_TYPE": "LYSO", "CHO": "0"}')

            lyso_dct = {'SMILES': 'O', 'OAP': 0, 'OCP': 0, 'DB': 0,
                        'OH': 0, 'KETO': 0, 'CHO': 0, 'COOH': 0, 'MOD_NUM': 0,
                        'FULL_SMILES': 'O', 'C_NUM': 0,
                        'FA_CHECKER': '0:0[0xDB,0xOH,0xKETO]<CHO@C0,COOH@C0>{OAP:0,OCP:0}',
                        'FA_ABBR': '0:0', 'FA_TYPE': 'LYSO', 'FA_JSON': lyso_json}

            _lyso_df = pd.DataFrame(lyso_dct, index=['0-lyso'])

            mod_sum_df = mod_sum_df.append(_lyso_df)

            return mod_sum_df

        return __bulk_oxidizer
    return _bulk_oxidizer


@bulk_oxidizer(TheoDB_Oxidizer)
def oxidizer(fa_link_dct):

    # fa_link_dct = {'FULL_smiles': usr_fa_smiles,
    #                'PRE_str': pre_ester_link_str,
    #                'POST_str': post_ester_link_str,
    #                'PRE_rgx': pre_ester_link_rgx_str,
    #                'POST_rgx': post_ester_link_rgx_str}

    db_info_dct = {}

    _usr_fa_smiles = fa_link_dct['FULL_smiles']
    _usr_fa_pre_str = fa_link_dct['PRE_str']
    _usr_fa_post_str = fa_link_dct['POST_str']
    _usr_fa_pre_rgx_str = fa_link_dct['PRE_rgx']
    _usr_fa_post_rgx_str = fa_link_dct['POST_rgx']
    _usr_fa_link_str = fa_link_dct['LINK_type']

    fa_rgx = re.compile(r'(%s)([C/\\=()]*)(%s)' % (_usr_fa_pre_rgx_str, _usr_fa_post_rgx_str))

    db_rgx = re.compile(r'(C/C[=]C\\)')

    # Construct regular expression for DB
    # if use FA: 'OC(CCCCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
    # get ('CCCCCC', 'C/C=C\\C/C=C\\C/C=C\\C/C=C\\', 'C/C=C\\', 'CCCC')
    fa_main_rgx = re.compile(r'(C*)((C[/]C[=]C\\){1,8})(C*)')

    fa_checker = re.match(fa_rgx, _usr_fa_smiles)
    if fa_checker:
        pre_fa_lst = fa_checker.groups()
        fa_main_str = pre_fa_lst[1]
        fa_main_checker = re.match(fa_main_rgx, fa_main_str)
        if fa_main_checker:
            pre_db_lst = fa_main_checker.groups()

            if len(pre_db_lst) >= 3:
                db_str = pre_db_lst[1]
                db_counter = db_str.count('C/C=C\\')
                # print('db_str', db_str, 'db_counter', db_counter)

                db_info_dct['DB_count'] = db_counter
                db_info_dct['DB_pre_part'] = ''.join([_usr_fa_pre_str, pre_db_lst[0]])
                db_info_dct['DB_main_part'] = pre_db_lst[1]
                db_info_dct['DB_post_part'] = ''.join([pre_db_lst[-1], _usr_fa_post_str])
                db_info_dct['DB_end_part'] = _usr_fa_post_str
                if (''.join([db_info_dct['DB_pre_part'], db_info_dct['DB_main_part'], db_info_dct['DB_post_part']])
                        == _usr_fa_smiles):
                    db_info_dct['DB_full_fa'] = _usr_fa_smiles
                    db_info_dct['DB_LINK_type'] = _usr_fa_link_str
                    db_info_dct['DB_C_count'] = _usr_fa_smiles.count('C')
                    db_info_dct['DB_start'] = db_info_dct['DB_pre_part'].count('C') + 2  # for that 2 C in DB pattern
                    db_info_dct['DB_omega'] = db_info_dct['DB_post_part'].count('C') + 1  # for that C in DB pattern

            else:
                pass
        else:
            db_info_dct['DB_end_part'] = _usr_fa_post_str
            db_info_dct['DB_full_fa'] = _usr_fa_smiles
            db_info_dct['DB_LINK_type'] = _usr_fa_link_str
            db_info_dct['DB_C_count'] = _usr_fa_smiles.count('C')
            db_info_dct['DB_count'] = 0
    else:
        db_info_dct['DB_end_part'] = _usr_fa_post_str
        db_info_dct['DB_full_fa'] = _usr_fa_smiles
        db_info_dct['DB_LINK_type'] = _usr_fa_link_str
        db_info_dct['DB_C_count'] = _usr_fa_smiles.count('C')
        db_info_dct['DB_count'] = 0

    return db_info_dct


def fa_link_filter(usr_fa_smiles):

    pre_ester_link_str = r'OC('
    pre_o_link_str = r'OCC'
    pre_p_link_str = r'O/C=C\C'

    post_ester_link_str = r'C)=O'
    post_o_link_str = r'C'
    post_p_link_str = r'C'

    pre_ester_link_rgx_str = r'OC[(]'
    pre_o_link_rgx_str = r'OCC'
    pre_p_link_rgx_str = r'O/C[=]C\\C'

    post_ester_link_rgx_str = r'C[)]=O'
    post_o_link_rgx_str = r'C'
    post_p_link_rgx_str = r'C'

    pre_ester_link_rgx = re.compile(pre_ester_link_rgx_str)
    pre_o_link_rgx = re.compile(pre_o_link_rgx_str)
    pre_p_link_rgx = re.compile(pre_p_link_rgx_str)

    # post_ester_link_rgx = re.compile(post_ester_link_rgx_str)
    # post_o_link_rgx = re.compile(post_o_link_rgx_str)
    # post_p_link_rgx = re.compile(post_p_link_rgx_str)

    ester_link_checker = re.match(pre_ester_link_rgx, usr_fa_smiles)
    o_link_checker = re.match(pre_o_link_rgx, usr_fa_smiles)
    p_link_checker = re.match(pre_p_link_rgx, usr_fa_smiles)
    
    if ester_link_checker:
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_ester_link_str, 'POST_str': post_ester_link_str,
                       'PRE_rgx': pre_ester_link_rgx_str, 'POST_rgx': post_ester_link_rgx_str, 'LINK_type': ''}
    elif o_link_checker:
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_o_link_str, 'POST_str': post_o_link_str,
                       'PRE_rgx': pre_o_link_rgx_str, 'POST_rgx': post_o_link_rgx_str, 'LINK_type': 'O-'}
    elif p_link_checker:
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_p_link_str, 'POST_str': post_p_link_str,
                       'PRE_rgx': pre_p_link_rgx_str, 'POST_rgx': post_p_link_rgx_str, 'LINK_type': 'P-'}
    else:
        fa_link_dct = {}

    # print(fa_link_dct)
    return fa_link_dct

# mod_table = 'ModConfig.csv'
# usr_fa = 'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
# # usr_fa = 'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
# # usr_fa = 'OCCCCCCCC/C=C\CCCCCCCC'
# # usr_fa = r'O/C=C\CCCCCC/C=C\CCCCCCCC'
# #
# fa_dct = fa_link_filter(usr_fa)
# mod_df = oxidizer(fa_dct, mod_table)
# print('mod_df', mod_df.shape)
