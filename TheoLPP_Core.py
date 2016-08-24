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

pl_class_use_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PS']

parser = PLParser()

pl_df = pd.read_excel(pl_table, sheetname=2)
fa_df = pd.read_csv(fa_table, index_col=0)

print(pl_df.head())
c_lst = []
for (_idx, _row) in pl_df.iterrows():

    _pl_abbr = _row['phospholipids']

    _pl_elem_lst = parser.get_composition(_pl_abbr)
    print ('PL composition ==>', _pl_elem_lst)
    _pl_hg_abbr = _pl_elem_lst[0]

    # get smiles from abbr

    if _pl_hg_abbr in pl_class_use_lst:
        c_lst.append(_pl_abbr)
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
            for (_sn2_idx, _sn2_row) in sn2_mod_sum_df.iterrows():
                _sn2_mod_smiles = _sn2_row['FULL_SMILES']

                if _pl_hg_abbr == 'PA':
                    _lpp_smiles = MergeBackLPP.pa_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr == 'PC':
                    _lpp_smiles = MergeBackLPP.pc_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr == 'PE':
                    _lpp_smiles = MergeBackLPP.pe_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr == 'PG':
                    _lpp_smiles = MergeBackLPP.pg_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr == 'PI':
                    _lpp_smiles = MergeBackLPP.pi_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr.upper() in ['PIP', 'PI4P']:
                    _lpp_smiles = MergeBackLPP.pi4p_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)
                elif _pl_hg_abbr == 'PS':
                    _lpp_smiles = MergeBackLPP.ps_lpp(_sn1_mod_smiles, _sn2_mod_smiles)
                    print(_pl_abbr,  _lpp_smiles)

print('==>==>==> %i of phospholipids processed==> ==> Finished !!!!' % len(c_lst))



    # fa_dct = fa_link_filter(usr_fa)
    # mod_df = oxidizer(fa_dct)
