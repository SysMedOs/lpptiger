# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova
# LPPtiger is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     LPPtiger repository: https://bitbucket.org/SysMedOs/lpptiger
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#

from __future__ import print_function
import re
import json

import numpy as np
import pandas as pd

from LibTheoLPP.AbbrGenerator import AbbrGenerator, fa_abbr_encode

from LibTheoLPP.TheoProstanes import IsoProstanOx
from LibTheoLPP.SMILESparser import SMILESparser


class TheoDB_Oxidizer:
    def __init__(self):
        print("Start to oxidize C=C -->")

    @staticmethod
    def mod_sum(usr_mod_table, usr_oxlevel):
        """

        :param str usr_mod_table:
        :return: the DataFrame of modification information
        :rtype: pd.DataFrame
        """
        mod_info_df = pd.read_csv(usr_mod_table, index_col=0, dtype={'OAP': np.int32})
        mod_info_df_t = mod_info_df.transpose()
        mod_info_df.loc['OXIDATIONLEVEL', :] = mod_info_df.loc['OXIDATIONLEVEL', :].astype(int)
        # print('mod_info_df_t')
        # print(mod_info_df_t)
        mod_info_df_t = mod_info_df_t[mod_info_df_t['OXIDATIONLEVEL'] < usr_oxlevel + 1]
        mod_info_df = mod_info_df_t.transpose()
        # print('mod_info_df')
        # print(mod_info_df)
        return mod_info_df


def oap_oxidizer(fa_dct, mod_info_df, ox_param_dct, mod_mode=0):
    max_mod = ox_param_dct['MAX_MOD']
    max_oh = ox_param_dct['MAX_MOD']
    max_keto = ox_param_dct['MAX_KETO']
    max_ooh = ox_param_dct['MAX_OOH']
    max_epoxy = ox_param_dct['MAX_EPOXY']

    db_count = fa_dct['DB_count']
    db_start_count = fa_dct['DB_start']

    db_count_i_lst = range(1, db_count + 1)
    db_count_s_lst = [str(db_x) for db_x in db_count_i_lst]
    unmod_lst = ['C/C=C/'] * db_count

    _unmod_dct = {}
    for _db_s in db_count_s_lst:
        _unmod_dct[_db_s] = 'C/C=C/'
        # _unmod_dct['MOD_' + _db_s] = ''.join([_db_s, '-', 'no_oxidation', '|'])
        _unmod_dct['MOD_' + _db_s] = ''
    _unmod_dct['MOD_COUNT'] = 0
    _unmod_dct['MOD_LIST'] = ''
    _mod_df = pd.DataFrame(_unmod_dct, index=[0])

    # get only the number of modifications
    if mod_mode == 0:

        oap_info_t_df = mod_info_df.transpose()
        oap_info_t_df = oap_info_t_df.reset_index()
        oap_info_t_df = oap_info_t_df.query('OAP > 0')
        print('oap_info_t_df', oap_info_t_df)
        mod_dct = {'OH': [], 'KETO': [], 'OOH': [], 'EPOXY': []}
        mod_idx_dct = {}
        for _idx, _mod_r in oap_info_t_df.iterrows():
            # after .t and reset_index, the mod name became column named 'index'
            mod_idx_dct[_mod_r['index']] = _idx
            for _k in mod_dct.keys():
                if _mod_r[_k] > 0:
                    _tmp_lst = mod_dct[_k]
                    _tmp_lst.append(_idx)
                    mod_dct[_k] = _tmp_lst

        oap_info_t_df['MAX_NUM'] = 0
        for _idx, _oap_row in oap_info_t_df.iterrows():

            _tmp_mod_count = _oap_row['MAX_NUM']
            if _oap_row['OH'] == 1:
                oap_info_t_df.set_value(_idx, 'MAX_NUM', _tmp_mod_count + max_mod)
                print('OH', oap_info_t_df)
            if _oap_row['KETO'] == 1:
                oap_info_t_df.set_value(_idx, 'MAX_NUM', _tmp_mod_count + max_keto)
            if _oap_row['OOH'] == 1:
                oap_info_t_df.set_value(_idx, 'MAX_NUM', _tmp_mod_count + max_ooh)
            if _oap_row['EPOXY'] == 1:
                oap_info_t_df.set_value(_idx, 'MAX_NUM', _tmp_mod_count + max_epoxy)
        print('oap_info_t_df')
        oap_info_t_df = oap_info_t_df.query('MAX_NUM > 0')
        print(oap_info_t_df)

        _tmp_oap_info_t_df = oap_info_t_df.copy()
        _tmp_db_mod_lst = []
        for db_i in db_count_i_lst:
            print('Now working on db number:', db_i)
            print(_tmp_oap_info_t_df.shape)
            _tmp_num_mod = _tmp_oap_info_t_df.shape[1]
            db_str = str(db_i)
            db_idx = db_i - 1

            # check if there is any mod left
            if _tmp_num_mod > 0:

                _rest_mod_df = _mod_df.query('MOD_COUNT == %i' % db_idx)

                for _mod_idx, mod_r in _rest_mod_df.iterrows():
                    if isinstance(mod_r['MOD_LIST'], str):
                        mod_lst = mod_r['MOD_LIST'].split('|')
                        print('mod_lst', mod_lst)
                        if len(mod_lst) == 0:
                            pass
                        else:
                            for _used_mod in mod_lst:
                                if len(_used_mod) > 0:
                                    pr_mod = _used_mod.split('-')
                                    _used_mod_name = pr_mod[1]
                                    print('_used_mod_name', _used_mod_name)
                                    _used_idx = mod_idx_dct[_used_mod_name]
                                    print('_used_idx', _used_idx)
                                    _tmp_max = _tmp_oap_info_t_df.get_value(_used_idx, 'MAX_NUM')
                                    _tmp_oap_info_t_df.set_value(_used_idx, 'MAX_NUM', _tmp_max - 1)
                            _tmp_oap_info_t_df = _tmp_oap_info_t_df.query('MAX_NUM > 0')
                    _tmp_mod_lst = _tmp_oap_info_t_df['SMILES'].tolist()
                    _tmp_mod_typ_lst = _tmp_oap_info_t_df['index'].tolist()
                    print(_tmp_mod_lst)
                    _tmp_mod_df = pd.DataFrame()
                    _tmp_mod_df[db_str] = _tmp_mod_lst
                    _tmp_db_mod = 'MOD_' + db_str
                    _tmp_db_mod_lst.append(_tmp_db_mod)
                    _tmp_mod_df[_tmp_db_mod] = _tmp_mod_typ_lst
                    _tmp_mod_df[_tmp_db_mod] = db_str + '-' + _tmp_mod_df[_tmp_db_mod] + '|'

                    # populate all possible mod
                    for db_s in db_count_s_lst:
                        if db_s == db_str:
                            pass
                        else:
                            # put rest part of db unchanged
                            _tmp_mod_df[db_s] = mod_r[db_s]
                    _tmp_mod_df['MOD_COUNT'] = db_i
                    _tmp_mod_df['MOD_LIST'] = _rest_mod_df['MOD_LIST'] + _tmp_mod_df[_tmp_db_mod]
                    _mod_df = _mod_df.append(_tmp_mod_df)
                    print(_mod_df)

                    # update the remaining number of oxidation sites
                    # _tmp_oap_info_t_df['MAX_NUM'] -= 1
                    # _tmp_oap_info_t_df = _tmp_oap_info_t_df.query('MAX_NUM > 0')
                    # print(_tmp_oap_info_t_df.shape)

    # get the exact position of each modification
    else:
        oap_info_t_df = pd.DataFrame()

    _mod_df = _mod_df.reset_index(drop=True)
    mod_df = fa_abbr_encode(fa_dct, oap_info_t_df, _mod_df, mode=0)

    print(mod_df)
    print(mod_df)

    return mod_df


# construct a decorator
def bulk_oxidizer(theodb_oxidizer_cls):
    def _bulk_oxidizer(ox_func):
        def __bulk_oxidizer(usr_fa_dct, usr_mod_table, isop_cfg, isopabbr_cfg,
                            oxlevel, ox_param_dct, prostane_mode, prostane_ox_mode):

            """

            :param usr_fa_dct:
            :param usr_mod_table:
            :param isop_cfg:
            :param isopabbr_cfg:
            :param oxlevel:
            :param ox_param_dct: ox_param_dct = {'MAX_MOD': int, 'MAX_KETO': int, 'MAX_OOH': int, 'MAX_EPOXY': int}
            :param prostane_mode:
            :param prostane_ox_mode:
            :return:
            """

            smi2formula = SMILESparser()

            max_mod = ox_param_dct['MAX_MOD']
            max_oh = ox_param_dct['MAX_MOD']
            max_keto = ox_param_dct['MAX_KETO']
            max_ooh = ox_param_dct['MAX_OOH']
            max_epoxy = ox_param_dct['MAX_EPOXY']

            nam_fields_last = ['OAP', 'OCP', 'DB', 'OH', 'KETO', 'CHO', 'COOH', 'OOH', 'EPOXY']

            fa_dct = ox_func(usr_fa_dct)
            # print('fa_dct', fa_dct)
            db_count = fa_dct['DB_count']
            fa_c_count = fa_dct['DB_C_count']
            # O- & P- link sn, the SMILES are different
            if fa_dct['DB_LINK_type'] in ['O-', 'P-']:
                ocp_end_part = ''
            else:
                ocp_end_part = ')=O'

            mod_info_df = theodb_oxidizer_cls.mod_sum(usr_mod_table, oxlevel)
            mod_typ_lst = mod_info_df.columns.tolist()
            mod_sum_df = pd.DataFrame()
            _def_frags_dct = {}
            for _mod_key in mod_typ_lst:

                if _mod_key not in ['aldehyde', 'aldehyde_short', 'carboxylic_acid', 'carboxylic_acid_short']:
                    _def_frags_dct[_mod_key] = ''.join([fa_dct['DB_pre_part'] +
                                                        mod_info_df.loc['FRAG', _mod_key] +
                                                        ocp_end_part])
                else:
                    _def_frags_dct[_mod_key] = ''

            _def_frags_df = pd.DataFrame(_def_frags_dct, index=['FRAG_SMILES'])
            mod_info_df = mod_info_df.append(_def_frags_df)

            del _def_frags_df, _def_frags_dct
            # print(mod_info_df.index.tolist())

            # change the data type from str to int
            mod_info_df.loc['OAP', :] = mod_info_df.loc['OAP', :].astype(int)
            mod_info_df.loc['OCP', :] = mod_info_df.loc['OCP', :].astype(int)
            mod_info_df.loc['DB', :] = mod_info_df.loc['DB', :].astype(int)
            mod_info_df.loc['OH', :] = mod_info_df.loc['OH', :].astype(int)
            mod_info_df.loc['KETO', :] = mod_info_df.loc['KETO', :].astype(int)
            mod_info_df.loc['CHO', :] = mod_info_df.loc['CHO', :].astype(int)
            mod_info_df.loc['COOH', :] = mod_info_df.loc['COOH', :].astype(int)
            mod_info_df.loc['OOH', :] = mod_info_df.loc['OOH', :].astype(int)
            mod_info_df.loc['EPOXY', :] = mod_info_df.loc['EPOXY', :].astype(int)
            # mod_info_df['FRAG_SMILES'] = ''

            if db_count > 0:

                # mod_df = oap_oxidizer(fa_dct, mod_info_df, ox_param_dct, mod_mode=0)
                # print(mod_df)

                db_range_lst = range(1, db_count + 1)
                # start oxidation
                for db_i in db_range_lst:

                    for _mod in mod_typ_lst:
                        _tmp_mod_lst = mod_sum_df.columns.tolist()
                        if len(_tmp_mod_lst) <= len(mod_typ_lst) and db_i == 1:
                            _mod_one = ''.join([str(db_i), '-', _mod])
                            mod_sum_df[_mod_one] = mod_info_df[_mod]
                            # add one more cleavage site for the first -OH bond
                            # if mod_sum_df.loc['FRAG_SMILES', _mod_one][-7:] == 'C(O))=O':
                            #     _mod_one_lst = ['["', mod_sum_df.loc['FRAG_SMILES', _mod_one][:-7], ocp_end_part,
                            #                     '","', mod_sum_df.loc['FRAG_SMILES', _mod_one], '"]']
                            #     _mod_one_json = ''.join(_mod_one_lst)
                            #     mod_sum_df.loc['FRAG_SMILES', _mod_one] = _mod_one_json
                            # elif mod_sum_df.loc['FRAG_SMILES', _mod_one][-4:] == 'C(O)':
                            #     _mod_one_lst = ['["', mod_sum_df.loc['FRAG_SMILES', _mod_one][:-4], '","',
                            #                     mod_sum_df.loc['FRAG_SMILES', _mod_one], '"]']
                            #     _mod_one_json = ''.join(_mod_one_lst)
                            #     mod_sum_df.loc['FRAG_SMILES', _mod_one] = _mod_one_json
                            # else:
                            #     _mod_one_json = ''.join(['["', mod_sum_df.loc['FRAG_SMILES', _mod_one], '"]'])
                            #     mod_sum_df.loc['FRAG_SMILES', _mod_one] = _mod_one_json

                        else:
                            for _tmp_mod in _tmp_mod_lst:
                                _tmp_mod = str(_tmp_mod)
                                __tmp_mod_lst = _tmp_mod.split('-')
                                # accelerate the speed by processing first few C=C for MODs
                                if int(__tmp_mod_lst[-2]) == db_i - 1:
                                    _ocp_checker = mod_sum_df.loc['OCP', _tmp_mod]
                                    if _ocp_checker == 0:
                                        _new_mod = ''.join([_tmp_mod, '-', str(db_i), '-', _mod])
                                        # this order is important for the SMILES generated
                                        # mod_sum_df.loc[:, _new_mod] = mod_sum_df[_tmp_mod] + mod_info_df[_mod]
                                        _tmp_mod_df = pd.DataFrame(mod_sum_df[_tmp_mod] + mod_info_df[_mod],
                                                                   columns=[_new_mod])

                                        # _tmp_mod_smiles = (fa_dct['DB_pre_part'] +
                                        #                    mod_sum_df.loc['SMILES', _tmp_mod] +
                                        #                    mod_info_df.loc['FRAG', _mod] + ocp_end_part)
                                        #
                                        # _tmp_frags_lst = json.loads(mod_sum_df.loc['FRAG_SMILES', _tmp_mod])

                                        # # Add one more cleavage site for OH
                                        # if _tmp_mod_smiles[-7:] == 'C(O))=O':
                                        #     _tmp_frags_lst.append(''.join([_tmp_mod_smiles[:-7], ocp_end_part]))
                                        # elif _tmp_mod_smiles[-4:] == 'C(O)':
                                        #     _tmp_frags_lst.append(_tmp_mod_smiles[:-4])

                                        # # Filter out full length OCPs
                                        # if db_i == db_count:
                                        #     if _tmp_mod_smiles[-4:] == 'C)=O':
                                        #         _tmp_frags_lst.append(_tmp_mod_smiles)
                                        #     elif _tmp_mod_smiles[-1:] == 'C':
                                        #         _tmp_frags_lst.append(_tmp_mod_smiles)
                                        # else:
                                        #     _tmp_frags_lst.append(_tmp_mod_smiles)
                                        #
                                        # _tmp_mod_df.loc['FRAG_SMILES', _new_mod] = json.dumps(_tmp_frags_lst)

                                        mod_sum_df.loc[:, _new_mod] = _tmp_mod_df

                                    else:
                                        # print('OCP', _ocp_checker, _tmp_mod)
                                        pass
                                else:
                                    pass
                mod_sum_df = mod_sum_df.transpose()
                # no more than one keto & max reduce 1 DB
                # _min_DB = db_count - 1
                # modification type and number control
                # mod_sum_df = mod_sum_df[(mod_sum_df.KETO < max_keto) & (mod_sum_df.DB >= _min_DB)]
                mod_sum_df['MOD_NUM'] = mod_sum_df['OAP'] + mod_sum_df['OCP']

                # filter the OCP and OAP. OAP should be full length
                mod_ocp_sum_df = mod_sum_df.query('OCP == 1')
                mod_ocp_sum_df = mod_ocp_sum_df.query('0 < MOD_NUM <= %d' % max_mod)
                # OAP should have all DB, thus MOD_NUM <= db_count
                mod_oap_sum_df = mod_sum_df.query('OCP == 0 and MOD_NUM <= %d' % db_count)
                mod_oap_sum_df = mod_oap_sum_df.query('0 < MOD_NUM <= %d' % max_mod)

                # the end of smiles is different for OCP
                mod_ocp_sum_idx_lst = mod_ocp_sum_df.index.tolist()
                # mod_ocp_sum_df.loc[mod_ocp_sum_idx_lst, 'FULL_SMILES'] = (fa_dct['DB_pre_part'] +
                #                                                           mod_ocp_sum_df['SMILES'] + ocp_end_part)
                # mod_ocp_sum_df.is_copy = False
                for _i_ocp, _r_ocp in mod_ocp_sum_df.iterrows():
                    mod_ocp_sum_df = mod_ocp_sum_df.set_value(_i_ocp, 'FULL_SMILES',
                                                              fa_dct['DB_pre_part'] + _r_ocp['SMILES']
                                                              + ocp_end_part)

                mod_oap_sum_idx_lst = mod_oap_sum_df.index.tolist()
                # mod_oap_sum_df.loc[mod_oap_sum_idx_lst, 'FULL_SMILES'] = (fa_dct['DB_pre_part'] +
                #                                                           mod_oap_sum_df['SMILES'] +
                #                                                           fa_dct['DB_post_part'])
                # mod_oap_sum_df.is_copy = False
                for _i_oap, _r_oap in mod_oap_sum_df.iterrows():
                    mod_oap_sum_df = mod_oap_sum_df.set_value(_i_oap, 'FULL_SMILES',
                                                              fa_dct['DB_pre_part'] + _r_oap['SMILES']
                                                              + fa_dct['DB_post_part'])

                mod_sum_df = mod_ocp_sum_df.append(mod_oap_sum_df)

                # select the LPP according to user settings of max modification types
                mod_max_ctrl_code = 'KETO <= %i and OOH<= %i' % (max_keto, max_ooh)
                mod_sum_df = mod_sum_df.query(mod_max_ctrl_code)

                # mod_sum_df = mod_sum_df.query('EPOXY <= %i' % max_epoxy)
                _full_smiles_lst = mod_sum_df['FULL_SMILES'].tolist()
                _c_num_lst = []
                for _smiles in _full_smiles_lst:
                    _c_num_lst.append(_smiles.count('C'))
                mod_sum_df.loc[:, 'C_NUM'] = _c_num_lst
                mod_sum_df.is_copy = False

                # OAP should stay the same carbon number
                _n_ocp_mod_sum_df = mod_sum_df.query('OCP == 1')
                _n_oap_mod_sum_df = mod_sum_df.query('OAP >= 1 and C_NUM == %d' % fa_c_count)
                mod_sum_df = _n_ocp_mod_sum_df.append(_n_oap_mod_sum_df)

                # mod_sum_df.to_csv('oxDB.csv')
                nam_fields_last.append('C_NUM')
                mod_sum_df[nam_fields_last] = mod_sum_df[nam_fields_last].astype(str)

                # here CHO@C & COOH@C is 1 or 0 for T/F
                mod_sum_df.loc[:, 'FA_CHECKER'] = (fa_dct['DB_LINK_type'] + mod_sum_df['C_NUM'] + ':'
                                                   + mod_sum_df['DB'] + '[' + mod_sum_df['DB'] + 'xDB,'
                                                   + mod_sum_df['OH'] + 'xOH,'
                                                   + mod_sum_df['KETO'] + 'xKETO,'
                                                   + mod_sum_df['OOH'] + 'xOOH,'
                                                   + mod_sum_df['EPOXY'] + 'xEPOXY]<CHO@C'
                                                   + mod_sum_df['CHO'] + ',COOH@C'
                                                   + mod_sum_df['COOH'] + '>{OAP:'
                                                   + mod_sum_df['OAP'] + ',OCP:'
                                                   + mod_sum_df['OCP'] + '}')
                abbr_gen = AbbrGenerator()
                mod_sum_df['FA_ABBR'] = ''
                mod_sum_df['FA_TYPE'] = ''
                mod_sum_df['FA_JSON'] = ''
                mod_sum_df['FA_FORMULA'] = ''

                for (_fa_idx, _fa_row) in mod_sum_df.iterrows():
                    _fa_code = str(_fa_row['FA_CHECKER'])
                    _fa_abbr, _fa_typ = abbr_gen.decode(_fa_code)
                    _fa_row['FA_ABBR'] = _fa_abbr
                    # print('_fa_code', _fa_code)
                    # print(_fa_typ, ' | ', _fa_abbr)
                    _fa_row['FA_TYPE'] = _fa_typ
                    _fa_row['FA_FORMULA'] = smi2formula.smiles2formula(_fa_row['FULL_SMILES'], charge='M')

                    # force all numbers to int, important for json coding!
                    _fa_row['FA_JSON'] = json.dumps({'LINK_TYPE': fa_dct['DB_LINK_type'],
                                                     'C': int(_fa_row['C_NUM']),
                                                     'DB': int(_fa_row['DB']),
                                                     'OH': int(_fa_row['OH']),
                                                     'KETO': int(_fa_row['KETO']),
                                                     'OOH': int(_fa_row['OOH']),
                                                     'EPOXY': int(_fa_row['EPOXY']),
                                                     'CHO': int(_fa_row['CHO']),
                                                     'COOH': int(_fa_row['COOH']),
                                                     'OAP': int(_fa_row['OAP']),
                                                     'OCP': int(_fa_row['OCP'])})

                if db_count >= 3 and prostane_mode == 1:
                    # isop_cfg = r'D:\theolpp\LibTheoLPP\IsoP_ModConfig.csv'
                    ox_prostane = IsoProstanOx(fa_dct, isop_cfg, isopabbr_cfg)

                    _prostane_lpp_dct = ox_prostane.get_isop_lpp()
                    _prostane_df = pd.DataFrame(_prostane_lpp_dct)
                    _prostane_df = _prostane_df.transpose()
                    for _idx_p, _prostane in _prostane_df.iterrows():
                        _prostane_df = _prostane_df.set_value(_idx_p, 'FA_FORMULA',
                                                              smi2formula.smiles2formula(_prostane['FULL_SMILES'],
                                                                                         charge='M'))

                    mod_sum_df = mod_sum_df.append(_prostane_df)
                    # print('mod_sum_df', mod_sum_df.shape)

            # print('mod_sum_df')
            # print(mod_sum_df)

            if fa_dct['DB_LINK_type'] == 'O-':
                _unmod_fa_abbr = 'O-%i:%i' % (fa_dct['DB_C_count'], db_count)
            elif fa_dct['DB_LINK_type'] == 'P-':
                _unmod_fa_abbr = 'P-%i:%i' % (fa_dct['DB_C_count'], db_count)
            elif fa_dct['DB_LINK_type'] == '':
                _unmod_fa_abbr = '%i:%i' % (fa_dct['DB_C_count'], db_count)
            else:
                _unmod_fa_abbr = '%i:%i' % (fa_dct['DB_C_count'], db_count)

            unmod_json = ('{"C": %i, "DB": %i, "CHO": 0, "EPOXY": 0, "OAP": 0, "OCP": 0, "COOH": 0, '
                          '"KETO": 0, "OH": 0, "OOH": 0, "LINK_TYPE": "%s"}'
                          % (fa_dct['DB_C_count'], db_count, fa_dct['DB_LINK_type']))

            unmod_dct = {'SMILES': fa_dct['DB_full_fa'], 'OAP': 0, 'OCP': 0, 'DB': db_count,
                         'OH': 0, 'KETO': 0, 'OOH': 0, 'EPOXY': 0, 'CHO': 0, 'COOH': 0, 'MOD_NUM': 0,
                         'FULL_SMILES': fa_dct['DB_full_fa'], 'C_NUM': fa_dct['DB_C_count'],
                         'FA_CHECKER': ('%i:%i[%ixDB,0xOH,0xKETO,0xOOH,0xEPOXY]'
                                        '<CHO@C0,COOH@C0>{OAP:0,OCP:0}'
                                        % (fa_dct['DB_C_count'], db_count, db_count)),
                         'FA_ABBR': _unmod_fa_abbr, 'FA_TYPE': 'UNMOD', 'FA_JSON': unmod_json, 'FRAG_SMILES': '[""]',
                         'FA_FORMULA': smi2formula.smiles2formula(fa_dct['DB_full_fa'], charge='M')}

            unmod_df = pd.DataFrame(unmod_dct, index=['0-no_oxidation'])

            lyso_json = ('{"C": 0, "DB": 0, "CHO": 0, "EPOXY": 0, "OAP": 0, "OCP": 0, "COOH": 0, '
                         '"KETO": 0, "OH": 0, "OOH": 0, "LINK_TYPE": "LYSO"}')

            lyso_dct = {'SMILES': 'O', 'OAP': 0, 'OCP': 0, 'DB': 0,
                        'OH': 0, 'KETO': 0, 'OOH': 0, 'EPOXY': 0, 'CHO': 0, 'COOH': 0, 'MOD_NUM': 0,
                        'FULL_SMILES': 'O', 'C_NUM': 0, 'FA_FORMULA': 'H2O',
                        'FA_CHECKER': '0:0[0xDB,0xOH,0xKETO,0xOOH,0xEPOXY]<CHO@C0,COOH@C0>{OAP:0,OCP:0}',
                        'FA_ABBR': '0:0', 'FA_TYPE': 'LYSO', 'FA_JSON': lyso_json, 'FRAG_SMILES': '[""]'}

            _lyso_df = pd.DataFrame(lyso_dct, index=['0-lyso'])

            mod_sum_df = mod_sum_df.append(unmod_df)
            mod_sum_df = mod_sum_df.append(_lyso_df)

            return mod_sum_df

        return __bulk_oxidizer

    return _bulk_oxidizer


@bulk_oxidizer(TheoDB_Oxidizer)
def oxidizer(fa_link_dct):
    """

    :param dict fa_link_dct: output from def fa_link_filter
    :return:
    """

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

    # db_rgx = re.compile(r'(C/C[=]C\\)')

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
            db_info_dct['DB_pre_part'] = fa_link_dct['FULL_smiles']
            db_info_dct['DB_end_part'] = _usr_fa_post_str
            db_info_dct['DB_full_fa'] = _usr_fa_smiles
            db_info_dct['DB_LINK_type'] = _usr_fa_link_str
            db_info_dct['DB_C_count'] = _usr_fa_smiles.count('C')
            db_info_dct['DB_count'] = 0
    else:
        db_info_dct['DB_pre_part'] = fa_link_dct['FULL_smiles']
        db_info_dct['DB_end_part'] = _usr_fa_post_str
        db_info_dct['DB_full_fa'] = _usr_fa_smiles
        db_info_dct['DB_LINK_type'] = _usr_fa_link_str
        db_info_dct['DB_C_count'] = _usr_fa_smiles.count('C')
        db_info_dct['DB_count'] = 0

    return db_info_dct


def fa_link_filter(usr_fa_smiles):
    """
    Find the P- / O- connections of sn1 / sn2 residues
    :param str usr_fa_smiles: the smiles code of sn1 / sn2 fatty acid
    :return dict fa_link_dct: the dictionary of P- / O- information
    """

    pre_ester_link_str = r'OC('
    pre_o_link_str = r'OCC'
    pre_p_link_str = r'O/C=C\C'

    post_ester_link_str = r'C)=O'
    post_o_link_str = r'C'
    post_p_link_str = r'C'

    post_ester_link_cho_str = r'C=O)=O'
    post_o_link_cho_str = r'C=O'
    post_p_link_cho_str = r'C=O'

    post_ester_link_cooh_str = r'C(O)=O)=O'
    post_o_link_cooh_str = r'C(O)=O'
    post_p_link_cooh_str = r'C(O)=O'

    pre_ester_link_rgx_str = r'OC[(]'
    pre_o_link_rgx_str = r'OCC'
    pre_p_link_rgx_str = r'O/C[=]C\\C'

    post_ester_link_rgx_str = r'C[)]=O'
    post_o_link_rgx_str = r'C'
    post_p_link_rgx_str = r'C'

    post_ester_link_cho_rgx_str = r'C=O[)]=O'
    post_o_link_cho_rgx_str = r'C=O'
    post_p_link_cho_rgx_str = r'C=O'

    post_ester_link_cooh_rgx_str = r'C[(]O[)]=O[)]=O'
    post_o_link_cooh_rgx_str = r'C[(]O[)][=]O'
    post_p_link_cooh_rgx_str = r'C[(]O[)][=]O'

    pre_ester_link_rgx = re.compile(pre_ester_link_rgx_str)
    pre_o_link_rgx = re.compile(pre_o_link_rgx_str)
    pre_p_link_rgx = re.compile(pre_p_link_rgx_str)

    post_ester_link_rgx = re.compile(r'.*%s' % post_ester_link_rgx_str)
    post_o_link_rgx = re.compile(r'.*%s' % post_o_link_rgx_str)
    post_p_link_rgx = re.compile(r'.*%s' % post_p_link_rgx_str)

    post_ester_link_cho_rgx = re.compile(r'.*%s' % post_ester_link_cho_rgx_str)
    post_o_link_cho_rgx = re.compile(r'.*%s' % post_o_link_cho_rgx_str)
    post_p_link_cho_rgx = re.compile(r'.*%s' % post_p_link_cho_rgx_str)

    post_ester_link_cooh_rgx = re.compile(r'.*%s' % post_ester_link_cooh_rgx_str)
    post_o_link_cooh_rgx = re.compile(r'.*%s' % post_o_link_cooh_rgx_str)
    post_p_link_cooh_rgx = re.compile(r'.*%s' % post_p_link_cooh_rgx_str)

    ester_link_checker = re.match(pre_ester_link_rgx, usr_fa_smiles)
    o_link_checker = re.match(pre_o_link_rgx, usr_fa_smiles)
    p_link_checker = re.match(pre_p_link_rgx, usr_fa_smiles)

    post_ester_link_checker = re.match(post_ester_link_rgx, usr_fa_smiles)
    post_o_link_checker = re.match(post_o_link_rgx, usr_fa_smiles)
    post_p_link_checker = re.match(post_p_link_rgx, usr_fa_smiles)

    post_ester_link_cho_checker = re.match(post_ester_link_cho_rgx, usr_fa_smiles)
    post_o_link_cho_checker = re.match(post_o_link_cho_rgx, usr_fa_smiles)
    post_p_link_cho_checker = re.match(post_p_link_cho_rgx, usr_fa_smiles)

    post_ester_link_cooh_checker = re.match(post_ester_link_cooh_rgx, usr_fa_smiles)
    post_o_link_cooh_checker = re.match(post_o_link_cooh_rgx, usr_fa_smiles)
    post_p_link_cooh_checker = re.match(post_p_link_cooh_rgx, usr_fa_smiles)

    post_rgx = ''
    post_str = ''

    if ester_link_checker:
        if post_ester_link_checker:
            post_str = post_ester_link_str
            post_rgx = post_ester_link_rgx_str
        elif post_ester_link_cho_checker:
            post_str = post_ester_link_cho_str
            post_rgx = post_ester_link_cho_rgx_str
        elif post_ester_link_cooh_checker:
            post_str = post_ester_link_cooh_str
            post_rgx = post_ester_link_cooh_rgx_str
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_ester_link_str, 'POST_str': post_str,
                       'PRE_rgx': pre_ester_link_rgx_str, 'POST_rgx': post_rgx, 'LINK_type': ''}
    elif o_link_checker:
        if post_o_link_checker:
            post_str = post_o_link_str
            post_rgx = post_o_link_rgx_str
        elif post_o_link_cho_checker:
            post_str = post_o_link_cho_str
            post_rgx = post_o_link_cho_rgx_str
        elif post_o_link_cooh_checker:
            post_str = post_o_link_cooh_str
            post_rgx = post_o_link_cooh_rgx_str
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_o_link_str, 'POST_str': post_str,
                       'PRE_rgx': pre_o_link_rgx_str, 'POST_rgx': post_rgx, 'LINK_type': 'O-'}
    elif p_link_checker:
        if post_p_link_checker:
            post_str = post_p_link_str
            post_rgx = post_p_link_rgx_str
        elif post_p_link_cho_checker:
            post_str = post_p_link_cho_str
            post_rgx = post_p_link_cho_rgx_str
        elif post_p_link_cooh_checker:
            post_str = post_p_link_cooh_str
            post_rgx = post_p_link_cooh_rgx_str
        fa_link_dct = {'FULL_smiles': usr_fa_smiles, 'PRE_str': pre_p_link_str, 'POST_str': post_str,
                       'PRE_rgx': pre_p_link_rgx_str, 'POST_rgx': post_rgx, 'LINK_type': 'P-'}
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
