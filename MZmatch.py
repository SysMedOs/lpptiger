# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import pandas as pd

from lpplibs.ExactMassCalc import Elem2Mass


class MZmatcher(object):

    def __init__(self, usr_csv):
        self.usr_csv = usr_csv

    def match2lib(self, csvlib, ppm=None):

        mzcalc = Elem2Mass()

        if 0 < ppm < 500:
            pass
        else:
            ppm = 10

        csvlib_df = pd.read_csv(csvlib, index_col=0, header=0, names=['PL_Abbr', 'Elem', 'ExactMass', 'mz_H', 'mz_Na', 'mz_K', 'SMILES'])
        print csvlib_df.head()

        usrcsv_df = pd.read_csv(self.usr_csv, names=['Target_m/z', 'H', 'Na', 'K', 'H_info', 'Na_info', 'K_info'])
        print usrcsv_df.tail()
        usrcsv_count = usrcsv_df.shape[0]

        found_usr_mz_lst = []
        found_lib_mz_lst = []
        found_lib_CHG_lst = []
        found_lib_name_lst = []
        found_lib_formula_lst = []

        for idx in range(0, usrcsv_count):
            usr_mz = usrcsv_df.iloc[idx, 0]
            # calc the m/z range based on user ppm, limit to 5 decimals
            usr_mz_range = (round(usr_mz * (1-ppm*0.000001), 5), round(usr_mz * (1+ppm*0.000001), 5))
            # print usr_mz, usr_mz_range

            # start query [M+H]+
            _query_H = str(usr_mz_range[0]) + ' < mz_H < ' + str(usr_mz_range[1])
            _found_df_H = csvlib_df.query(_query_H)
            if _found_df_H.shape[0] > 0:
                # print usr_mz, _found_df_H
                usrcsv_df.loc[idx, 'H'] = 'H'
                info_H = ''
                _idx_lst_H = _found_df_H.index.tolist()

                # print _idx_lst, _found_df_H.loc[_idx_lst, 'PL_Abbr'], _found_df_H.loc[_idx_lst, 'mz_H']
                for _idx_H in _idx_lst_H:
                    found_usr_mz_lst.append(usr_mz)
                    found_lib_CHG_lst.append('H')
                    found_lib_mz_lst.append(str(_found_df_H.loc[_idx_H, 'mz_H']))
                    found_lib_name_lst.append(str(_found_df_H.loc[_idx_H, 'PL_Abbr']))
                    _formula_H_raw = str(_found_df_H.loc[_idx_H, 'Elem']) + 'H'
                    _formula_H = mzcalc.get_elem(_formula_H_raw)
                    found_lib_formula_lst.append(_formula_H)
                    if info_H == '':
                        info_H += (str(_found_df_H.loc[_idx_H, 'PL_Abbr']) + ': ' +
                                  str(_found_df_H.loc[_idx_H, 'mz_H']))
                    else:
                        info_H += ('| ' + str(_found_df_H.loc[_idx_H, 'PL_Abbr']) + ': ' +
                                   str(_found_df_H.loc[_idx_H, 'mz_H']))
                print info_H
                usrcsv_df.loc[idx, 'H_info'] = info_H
            # start query [M+Na]+
            _query_Na = str(usr_mz_range[0]) + ' < mz_Na < ' + str(usr_mz_range[1])
            _found_df_Na = csvlib_df.query(_query_Na)
            if _found_df_Na.shape[0] > 0:
                # print usr_mz, _found_df_Na
                usrcsv_df.loc[idx, 'Na'] = 'Na'
                info_Na = ''
                _idx_lst_Na = _found_df_Na.index.tolist()
                # print _idx_lst, _found_df_H.loc[_idx_lst, 'PL_Abbr'], _found_df_H.loc[_idx_lst, 'mz_Na']
                for _idx_Na in _idx_lst_Na:
                    found_usr_mz_lst.append(usr_mz)
                    found_lib_CHG_lst.append('Na')
                    found_lib_mz_lst.append(str(_found_df_Na.loc[_idx_Na, 'mz_Na']))
                    found_lib_name_lst.append(str(_found_df_Na.loc[_idx_Na, 'PL_Abbr']))
                    _formula_Na = str(_found_df_Na.loc[_idx_Na, 'Elem']) + 'Na'
                    found_lib_formula_lst.append(_formula_Na)
                    if info_Na == '':
                        info_Na += (str(_found_df_Na.loc[_idx_Na, 'PL_Abbr']) + ': ' +
                                    str(_found_df_Na.loc[_idx_Na, 'mz_Na']))
                    else:
                        info_Na += ('| ' + str(_found_df_Na.loc[_idx_Na, 'PL_Abbr']) + ': ' +
                                    str(_found_df_Na.loc[_idx_Na, 'mz_Na']))
                print info_Na
                usrcsv_df.loc[idx, 'Na_info'] = info_Na

            # start query [M+K]+
            _query_K = str(usr_mz_range[0]) + ' < mz_K < ' + str(usr_mz_range[1])
            _found_df_K = csvlib_df.query(_query_K)

            if _found_df_K.shape[0] > 0:
                # print _found_df_K
                # print usr_mz, _found_df_K
                usrcsv_df.loc[idx, 'K'] = 'K'
                info_K = ''
                _idx_lst_K = _found_df_K.index.tolist()
                # print _idx_lst, _found_df_H.loc[_idx_lst, 'PL_Abbr'], _found_df_H.loc[_idx_lst, 'mz_K']
                print _idx_lst_K
                for _idx_K in _idx_lst_K:
                    found_usr_mz_lst.append(usr_mz)
                    found_lib_CHG_lst.append('K')
                    found_lib_mz_lst.append(str(_found_df_K.loc[_idx_K, 'mz_K']))
                    found_lib_name_lst.append(str(_found_df_K.loc[_idx_K, 'PL_Abbr']))
                    _formula_K = str(_found_df_K.loc[_idx_K, 'Elem']) + 'K'
                    found_lib_formula_lst.append(_formula_K)
                    if info_K == '':
                        info_K += (str(_found_df_K.loc[_idx_K, 'PL_Abbr']) + ': ' +
                                   str(_found_df_K.loc[_idx_K, 'mz_K']))
                    else:
                        info_K += ('| ' + str(_found_df_K.loc[_idx_K, 'PL_Abbr']) + ': ' +
                                   str(_found_df_K.loc[_idx_K, 'mz_K']))
                print info_K
                usrcsv_df.loc[idx, 'K_info'] = info_K

        usrcsv_df.to_csv('user_marked_list.csv')

        found_dct = {'User_mz': found_usr_mz_lst, 'Found_mz': found_lib_mz_lst, 'Found_charge': found_lib_CHG_lst,
                     'Found_formula': found_lib_formula_lst, 'info': found_lib_name_lst}

        found_sum_df = pd.DataFrame(data=found_dct, columns=['User_mz', 'Found_mz', 'Found_charge',
                                                             'Found_formula', 'info'])
        found_sum_df.to_csv('Found_mz_table.csv')


matcher = MZmatcher('PotenzielleMarkerliste_MP19.csv')
matcher.match2lib('SDF_summary_table_lite.csv', ppm=2)

print 'Finished'
