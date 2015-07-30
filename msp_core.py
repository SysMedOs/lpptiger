# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import pandas as pd

from LPPfrag import TheoFrag

csv = r'PC(16-0_18-2(9-Z;12-Z)).csv'

df = pd.read_csv(csv, index_col=0)

idx_lst = range(len(df.index))

frag_obj = TheoFrag()
sum_dct = {}
for _idx in idx_lst:
    s = df.loc[_idx, 'SMILES']
    d = df.loc[_idx, 'PL_Abbr']
    _dct = frag_obj.smiles2frag(s, d, chargelist=['[M+H]+', '[M+Na]+', '[M-H2O+H]+', '[M-H2O+Na]+'], plclass='PC')

    for _key in _dct.keys():
        sum_dct[_key] = _dct[_key]

output = 'oxPLPC_frag.msp'
frag_obj.frag2msp(sum_dct, output)

