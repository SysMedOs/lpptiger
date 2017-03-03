# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPsmi.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

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
        try:
            _msp_json = _mol.GetProp('MSP_JSON')
        except KeyError:
            _msp_json = ''
        try:
            _fingerprint = _mol.GetProp('FINGERPRINT')
        except KeyError:
            _fingerprint = ''
        _full_smiles = _mol.GetProp('LPP_SMILES')
        _sn1_smiles = _mol.GetProp('SN1_SMILES')
        _sn2_smiles = _mol.GetProp('SN2_SMILES')

        _mol_info_dct = {'Abbreviation': _id, 'Class': _class, 'FORMULA_NEUTRAL': _formula,
                         'EXACT_MASS': _exactmass, 'MSP_JSON': _msp_json, 'FINGERPRINT': _fingerprint,
                         'LPP_SMILES': _full_smiles, 'SN1_SMILES': _sn1_smiles, 'SN2_SMILES': _sn2_smiles}
        for _charge in ['[M-H]-', '[M+HCOO]-']:
            if _charge in _pr_info_dct.keys():
                _mol_info_dct[_charge + '_FORMULA'] = _pr_info_dct[_charge][0]
                _mol_info_dct[_charge + '_MZ'] = _pr_info_dct[_charge][1]
            else:
                _mol_info_dct[_charge + '_FORMULA'] = np.NaN
                _mol_info_dct[_charge + '_MZ'] = np.NaN

        mol_info_dct[_id] = _mol_info_dct

    mol_info_df = pd.DataFrame(data=mol_info_dct)

    mol_info_df = mol_info_df.transpose()
    mol_info_df = mol_info_df[['Class', 'Abbreviation', 'FORMULA_NEUTRAL', 'EXACT_MASS',
                               '[M-H]-_FORMULA', '[M-H]-_MZ', '[M+HCOO]-_FORMULA', '[M+HCOO]-_MZ', 'MSP_JSON',
                               'FINGERPRINT', 'LPP_SMILES', 'SN1_SMILES', 'SN2_SMILES']]
    mol_info_df = mol_info_df.sort_values(by='EXACT_MASS')
    mol_info_df.to_excel(save_path, index=False)

    print('saved!')


def sdf2sum_fa(usr_sdf, save_path):

    mzcalc = MZcalc()

    usr_sdf = str(os.path.abspath(usr_sdf))

    mol_suppl = Chem.SDMolSupplier(usr_sdf)

    fa_info_lst = []

    for _mol in mol_suppl:

        _sn1_info_lst = [_mol.GetProp('SN1_ABBR'), _mol.GetProp('SN1_FORMULA'), _mol.GetProp('SN1_JSON')]
        _sn2_info_lst = [_mol.GetProp('SN2_ABBR'), _mol.GetProp('SN2_FORMULA'), _mol.GetProp('SN2_JSON')]

        if _sn1_info_lst in fa_info_lst:
            pass
        else:
            fa_info_lst.append(_sn1_info_lst)

        if _sn2_info_lst in fa_info_lst:
            pass
        else:
            fa_info_lst.append(_sn2_info_lst)

    mol_info_df = pd.DataFrame(data=fa_info_lst, columns=['FA', 'elem', 'INFO_JSON'])

    for i, r in mol_info_df.iterrows():
        _formula = r['elem']
        _info_dct = json.loads(r['INFO_JSON'])
        _link = _info_dct['LINK_TYPE']
        if _link in ['O', 'O-']:
            _link = 'O'
        elif _link in ['P', 'P-']:
            _link = 'P'
        else:
            _link = 'A'

        exact_mz = mzcalc.get_mono_mz(_formula, charge='')
        frag_mz = mzcalc.get_mono_mz(_formula, charge='[M-H]-')

        mol_info_df = mol_info_df.set_value(i, 'C', _info_dct['C'])
        mol_info_df = mol_info_df.set_value(i, 'DB', _info_dct['DB'])
        mol_info_df = mol_info_df.set_value(i, 'Link', _link)
        mol_info_df = mol_info_df.set_value(i, 'mass', exact_mz)
        mol_info_df = mol_info_df.set_value(i, '[M-H]-', frag_mz)
        mol_info_df = mol_info_df.set_value(i, '[M-H2O-H]-', exact_mz - 18.010565)
        mol_info_df = mol_info_df.set_value(i, 'NL-H2O', exact_mz - 18.010565)
    mol_info_df = mol_info_df.sort_values(by='mass')
    mol_info_df = mol_info_df[['FA', 'Link', 'C', 'DB', 'elem', 'mass', '[M-H]-', '[M-H2O-H]-', 'NL-H2O']]
    mol_info_df.to_excel(save_path)

    print('saved!')

if __name__ == '__main__':

    sdf = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.sdf'
    output = r'D:\theolpp\PX_18-1-2_max_1keto_1lessDB_FRAG.xlsx'

    sdf2xlsx(sdf, output)

