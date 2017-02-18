# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import print_function

import json
import os

import pandas as pd
import numpy as np

from rdkit import Chem

from ExactMassCalc import MZcalc


def sdf2xlsx(usr_sdf, save_path):

    usr_sdf = str(os.path.abspath(usr_sdf))

    mol_suppl = Chem.SDMolSupplier(usr_sdf)

    mol_info_dct = {}

    for _mol in mol_suppl:

        _id = _mol.GetProp('LM_ID')
        _class = _mol.GetProp('LPP_CLASS')
        _formula = _mol.GetProp('FORMULA')
        _exactmass = _mol.GetProp('EXACT_MASS')
        _pr_info_dct = json.loads(_mol.GetProp('PRECURSOR_JSON'))

        _mol_info_dct = {'CLASS': _class, 'FORMULA': _formula, 'EXACT_MASS': _exactmass}
        for _charge in ['[M-H]-', '[M+HCOO]-']:
            if _charge in _pr_info_dct.keys():
                _mol_info_dct[_charge] = _pr_info_dct[_charge][1]
            else:
                _mol_info_dct[_charge] = np.NaN

        mol_info_dct[_id] = _mol_info_dct

    mol_info_df = pd.DataFrame(data=mol_info_dct)

    mol_info_df = mol_info_df.transpose()
    mol_info_df = mol_info_df.sort_values(by='EXACT_MASS')
    mol_info_df.to_excel(save_path)

    print('saved!')


def sdf2sum_fa(usr_sdf, save_path):

    mzcalc = MZcalc()

    usr_sdf = str(os.path.abspath(usr_sdf))

    mol_suppl = Chem.SDMolSupplier(usr_sdf)

    fa_info_lst = []

    for _mol in mol_suppl:

        _sn1_info_lst = [_mol.GetProp('SN1_ABBR'), _mol.GetProp('SN1_FORMULA')]
        _sn2_info_lst = [_mol.GetProp('SN2_ABBR'), _mol.GetProp('SN2_FORMULA')]

        if _sn1_info_lst in fa_info_lst:
            pass
        else:
            fa_info_lst.append(_sn1_info_lst)

        if _sn2_info_lst in fa_info_lst:
            pass
        else:
            fa_info_lst.append(_sn2_info_lst)

    mol_info_df = pd.DataFrame(data=fa_info_lst, columns=['FA_ABBR', 'FORMULA'])

    for i, r in mol_info_df.iterrows():
        _formula = r['FORMULA']
        mol_info_df = mol_info_df.set_value(i, 'EXACT_MASS', mzcalc.get_mono_mz(_formula, charge=''))
        mol_info_df = mol_info_df.set_value(i, '[M-H]-', mzcalc.get_mono_mz(_formula, charge='[M-H]-'))
    mol_info_df = mol_info_df.sort_values(by='EXACT_MASS')
    mol_info_df.to_excel(save_path)

    print('saved!')

if __name__ == '__main__':

    sdf = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.sdf'
    output = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.xlsx'

    sdf2xlsx(sdf, output)

