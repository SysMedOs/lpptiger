# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
# import re

import pandas as pd

from lpplibs.PLParser import PLParser
from lpplibs.DBoxTheo import fa_link_filter, oxidizer
from lpplibs import MergeBackLPP


pl_table = './lpplibs/CM_NormalLipids.xlsx'
fa_table = './lpplibs/FA_list.csv'
mod_table = './lpplibs/ModConfig.csv'

pl_class_use_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']

parser = PLParser()

pl_df = pd.read_excel(pl_table, sheetname=2)
fa_df = pd.read_csv(fa_table, index_col=0)

print(pl_df.head())
c_lst = []

sum_theo_lpp_dct = {}
for (_idx, _row) in pl_df.iterrows():

    _pl_abbr = _row['phospholipids']

    _pl_elem_lst = parser.get_composition(_pl_abbr)
    print ('PL composition ==>', _pl_elem_lst)
    _pl_hg_abbr = _pl_elem_lst[0]

    # get smiles from abbr

    if _pl_hg_abbr in pl_class_use_lst:
        c_lst.append(_pl_abbr)

        # prepare output
        _pl_lpp_df = pd.DataFrame()

        print('Start oxidation of ==>', _pl_abbr)
        _pl_sn1_abbr = _pl_elem_lst[1]
        _pl_sn2_abbr = _pl_elem_lst[2]
        _pl_sn1_smiles = fa_df.loc[_pl_sn1_abbr, 'SMILES']
        _pl_sn2_smiles = fa_df.loc[_pl_sn2_abbr, 'SMILES']
        print('sn1 =>', _pl_sn1_smiles, '|| sn2 =>', _pl_sn1_smiles)
        sn1_link_dct = fa_link_filter(_pl_sn1_smiles)
        sn1_mod_sum_df = oxidizer(sn1_link_dct, mod_table)

        sn2_link_dct = fa_link_filter(_pl_sn2_smiles)
        sn2_mod_sum_df = oxidizer(sn2_link_dct, mod_table)

        for (_sn1_idx, _sn1_row) in sn1_mod_sum_df.iterrows():
            _sn1_mod_smiles = _sn1_row['FULL_SMILES']
            # print(_sn1_row)
            for (_sn2_idx, _sn2_row) in sn2_mod_sum_df.iterrows():
                _sn2_mod_smiles = _sn2_row['FULL_SMILES']
                # print(_sn2_row)

                # TODO(zhixu.ni@uni-leipzig.de): take more info of each mod from sn1 & sn2
                # _lpp_info_df = pd.DataFrame(data={'sn1': _sn1_row, 'sn2': _sn2_row})

                _lpp_smiles = MergeBackLPP.pl_lpp(_pl_hg_abbr, _sn1_mod_smiles, _sn2_mod_smiles)
                _lpp_info_dct = {'LPP_ORIGIN': _pl_abbr, 'LPP_SMILES': _lpp_smiles, 'LPP_CLASS': _pl_hg_abbr,
                            'SN1_SMILES': _sn1_mod_smiles, 'SN2_SMILES': _sn2_mod_smiles,
                            'SN1_INFO': _sn1_row['FA_CHECKER'], 'SN2_INFO': _sn2_row['FA_CHECKER']}
                _lpp_id_str = ''.join([_pl_hg_abbr, '(', _sn1_row['FA_CHECKER'], '/', _sn2_row['FA_CHECKER'], ')'])
                _lpp_info_se = pd.Series(data=_lpp_info_dct)
                _pl_lpp_df[_lpp_id_str] = _lpp_info_se

                del _lpp_info_dct, _lpp_info_se

        _pl_lpp_df = _pl_lpp_df.transpose()
        print('_pl_lpp_df', _pl_lpp_df.shape)
        print(_pl_lpp_df.head(5))
        sum_theo_lpp_dct[_pl_abbr] = _pl_lpp_df

sum_theo_lpp_pl = pd.Panel(data=sum_theo_lpp_dct)
print(sum_theo_lpp_pl.shape)

print('==>==>==> %i of phospholipids processed==> ==> Finished !!!!' % len(c_lst))
