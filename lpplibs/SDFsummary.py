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
        for _charge in ['[M-H]-', '[M+FA-H]-']:
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

if __name__ == '__main__':

    sdf = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.sdf'
    output = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.xlsx'

    sdf2xlsx(sdf, output)

