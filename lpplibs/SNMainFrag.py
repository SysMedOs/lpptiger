# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re
import json

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


class SNFrag(object):
    def __init__(self, pl_class, frag_score_list):
        self.pl_type = pl_class
        pl_frag_df = pd.read_excel(frag_score_list)
        pl_frag_df = pl_frag_df[pl_frag_df['PL_CLASS'] == self.pl_type]
        pl_frag_df['sn2_FA'] = pl_frag_df['sn2_FA'].astype('int')
        pl_frag_df['sn2_DB'] = pl_frag_df['sn2_DB'].astype('int')
        pl_frag_df['OH'] = pl_frag_df['OH'].astype('int')
        pl_frag_df['KETO'] = pl_frag_df['KETO'].astype('int')
        self.pl_frag_df = pl_frag_df

    def calc_frags(self, lpp_info_dct, mode='neg'):

        print(lpp_info_dct['LM_ID'])

        _lpp_type = lpp_info_dct['LPP_CLASS']
        _lpp_full_smi = lpp_info_dct['LPP_SMILES']
        _sn1_smi = lpp_info_dct['SN1_SMILES']
        _sn2_smi = lpp_info_dct['SN2_SMILES']
        _sn1_mod_dct = json.loads(lpp_info_dct['SN1_JSON'])
        _sn2_mod_dct = json.loads(lpp_info_dct['SN2_JSON'])

        df_header_lst = self.pl_frag_df.columns.tolist()
        frag_type_lst = []
        for _header in df_header_lst:
            if _header[-1] == '-':
                frag_type_lst.append(_header)

        # construct the return dict
        _frag_dct = {}

        # filter out both sn1 and sn2 modified species
        if _sn1_mod_dct['OAP'] + _sn1_mod_dct['OCP'] > 0 and _sn2_mod_dct['OAP'] + _sn2_mod_dct['OCP'] > 0:
            print('====>>>>>sn1 and sn2 both modified!! currently not supported!! skip now ====>>>>>')

        elif _sn1_mod_dct['LINK_TYPE'] in ['O-', 'P-'] or _sn2_mod_dct['LINK_TYPE'] in ['O-', 'P-']:
            print('====>>>>>sn1 or sn2 O-/P- link!! currently not supported!! skip now ====>>>>>')
        else:

            # get fragment scores in pandas series
            _score_se = self.get_score_type(_sn2_mod_dct, _lpp_type)
            # print(_score_se)

            if _lpp_type == 'PC':
                _frag_dct['[M+FA-H]-'] = self.calc_mz(_lpp_full_smi) + (int(_score_se.loc['[M+FA-H]-']),)
                _frag_dct['[M-H]-'] = self.calc_mz(_lpp_full_smi) + (int(_score_se.loc['[M-H]-']),)
            else:
                _frag_dct['[M-H]-'] = self.calc_mz(_lpp_full_smi) + (int(_score_se.loc['[M-H]-']),)

            if '[sn1-H]-' in frag_type_lst:
                _sn_info_tp = self.calc_mz(_sn1_smi)
                if _sn1_mod_dct['LINK_TYPE'] == 'LYSO':
                    pass
                else:
                    _frag_dct['[sn1-H]-'] = _sn_info_tp + (int(_score_se.loc['[sn1-H]-']),)

            if '[sn2-H]-' in frag_type_lst:
                _frag_dct['[sn2-H]-'] = self.calc_mz(_sn2_smi) + (int(_score_se.loc['[sn2-H]-']),)

                print(_frag_dct)

    def get_score_type(self, sn_mod_dct, lpp_class):

        # json dict format as below
        # '{"C": 16, "KETO": 0, "OH": 0, "OAP": 1, "OCP": 0, "COOH": 0, "DB": 1, "LINK_TYPE": "", "CHO": 0}'
        sn_c = sn_mod_dct['C']
        sn_db = sn_mod_dct['DB']
        sn_oh = sn_mod_dct['OH']
        sn_keto = sn_mod_dct['KETO']
        sn_oap = sn_mod_dct['OAP']
        sn_ocp = sn_mod_dct['OCP']
        sn_cho = sn_mod_dct['CHO']
        sn_cooh = sn_mod_dct['COOH']
        # sn_link = sn_mod_dct['LINK_TYPE']

        if sn_oap == 1:
            _lpp_score_df = self.pl_frag_df[(self.pl_frag_df['TYPE'] == 'OAP') &
                                            (self.pl_frag_df['PL_CLASS'] == lpp_class)]
            for (_lpp_idx, _lpp_r) in _lpp_score_df.iterrows():
                if sn_c == _lpp_r['sn2_FA'] and sn_db == _lpp_r['sn2_DB']:
                    if sn_oh == _lpp_r['OH'] and sn_keto == _lpp_r['KETO']:
                        return _lpp_r
                    else:
                        pass
                else:
                    pass
            # if nothing match, set -1 for any number of C and DB
            _lpp_r = _lpp_score_df[(_lpp_score_df['sn2_FA'] == -1) & (_lpp_score_df['sn2_DB'] == -1)]
            return _lpp_r.iloc[0, :]

        elif sn_ocp == 1 and sn_cho == 1:
            _lpp_score_df = self.pl_frag_df[self.pl_frag_df['TYPE'] == 'OCP_CHO']
            for (_lpp_idx, _lpp_r) in _lpp_score_df.iterrows():
                if sn_c == _lpp_r['sn2_FA'] and sn_db == _lpp_r['sn2_DB']:
                    if sn_oh == _lpp_r['OH'] and sn_keto == _lpp_r['KETO']:
                        return _lpp_r
                    else:
                        pass
                else:
                    pass
            _lpp_r = _lpp_score_df[(_lpp_score_df['sn2_FA'] == -1) & (_lpp_score_df['sn2_DB'] == -1)]
            return _lpp_r.iloc[0, :]

        elif sn_ocp == 1 and sn_cooh == 1:
            _lpp_score_df = self.pl_frag_df[self.pl_frag_df['TYPE'] == 'OCP_COOH']
            for (_lpp_idx, _lpp_r) in _lpp_score_df.iterrows():
                if sn_c == _lpp_r['sn2_FA'] and sn_db == _lpp_r['sn2_DB']:
                    if sn_oh == _lpp_r['OH'] and sn_keto == _lpp_r['KETO']:
                        return _lpp_r
                    else:
                        pass
                else:
                    pass
            _lpp_r = _lpp_score_df[(_lpp_score_df['sn2_FA'] == -1) & (_lpp_score_df['sn2_DB'] == -1)]
            return _lpp_r.iloc[0, :]

    def get_frag_scores(self):
        pass

    @staticmethod
    def calc_mz(smi, charge='[M-H]-'):

        charge_mz_dct = {'[M+H]+': 1.007825, '[M+Na]+': 22.989770, '[M+K]+': 38.963708, '[M+NH4]+': 18.034374,
                         '[M-H]-': -1.007825, '[M+FA-H]-': 44.997654, '[M+HCOO]-': 44.997654}

        charge_elem_dct = {'[M+H]+': {'H': 1}, '[M+Na]+': {'Na': 1}, '[M+K]+': {'Na': 1}, '[M+NH4]+': {'H': 4, 'N': 1},
                           '[M-H]-': {'H': -1}, '[M+FA-H]-': {'H': 1, 'C': 1, 'O': 2},
                           '[M+HCOO]-': {'H': 1, 'C': 1, 'O': 2}}

        if charge in charge_mz_dct.keys() and charge in charge_elem_dct.keys():
            pass
        else:
            charge = '[M-H]-'

        _mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(_mol)

        _formula = rdMolDescriptors.CalcMolFormula(_mol)
        _exactmass = rdMolDescriptors.CalcExactMolWt(_mol)

        _elem_rgx = re.compile(r'[A-Z][a-z]?[0-9]{0,3}')
        _elem_bulk_lst = re.findall(_elem_rgx, _formula)
        _elem_num_rgx = re.compile(r'([A-Z][a-z]?)([0-9]{0,3})')
        _elem_dct = {}
        for _elem in _elem_bulk_lst:
            _elem_checker = re.match(_elem_num_rgx, _elem)
            _elem_lst = _elem_checker.groups()
            try:
                _elem_dct[_elem_lst[0]] = int(_elem_lst[1])
            # e.g. only 1 N or P will get '' for _elem_lst[1]
            except ValueError:
                _elem_dct[_elem_lst[0]] = 1

        # print(_elem_bulk_lst, _exactmass)
        _elem_charged_dct = {}

        # get sum keys form both dict
        _charged_keys_lst = set(sum([_elem_dct.keys(), charge_elem_dct[charge].keys()], []))
        for _key in _charged_keys_lst:
            _elem_charged_dct[_key] = _elem_dct.get(_key, 0) + charge_elem_dct[charge].get(_key, 0)

        charge_mz = _exactmass + charge_mz_dct[charge]
        elem_order_lst = ['C', 'H', 'N', 'O', 'P', 'S', 'Na', 'K']
        _charged_elem = ''
        for _e in elem_order_lst:
            if _e in _charged_keys_lst:
                _charged_elem += _e
                if _elem_charged_dct[_e] > 1:
                    _charged_elem += str(_elem_charged_dct[_e])
                else:
                    pass
        if charge in ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+']:
            _charged_elem += '+'
        elif charge in ['[M-H]-', '[M+FA-H]-', '[M+HCOO]-']:
            _charged_elem += '-'

        charged_info = (_charged_elem, round(charge_mz, 4))

        return charged_info


# usr_smi = r'OC(CCCCC(O)/C=C/C(O)/C=C/C(O)/C=C/CC(O)=O)=O'
usr_smi = r'OP(OCCN)(OCC([H])(OC(CCCCCCC/C=C/C/C=C/CCCCC)=O)COC(CCCCCCC(O)CCC/C=C/CCCCC)=O)=O'

lpp_dct1 = {
    'LPP_FRAG': '["OP(OCCN)(OCC([H])(OC(CCCC)=O)CO/C=C\\\\CCCCCCCCCCCCCC)=O", "OP(OCCN)(OCC([H])(OC(CCCC/C=C/C)=O)CO/C=C\\\\CCCCCCCCCCCCCC)=O", "OP(OCCN)(OCC([H])(OC(CCCC/C=C/C/C=C/)=O)CO/C=C\\\\CCCCCCCCCCCCCC)=O", "OP(OCCN)(OCC([H])(OC(CCCC/C=C/C/C=C/C(O))=O)CO/C=C\\\\CCCCCCCCCCCCCC)=O", "OP(OCCN)(OCC([H])(OC(CCCC/C=C/C/C=C/C(O)/C=C/CC(O)=O)=O)CO/C=C\\\\CCCCCCCCCCCCCC)=O"]',
    'SN2_FRAGS': '["OC(CCCC)=O", "OC(CCCC/C=C/C)=O", "OC(CCCC/C=C/C/C=C/)=O", "OC(CCCC/C=C/C/C=C/C(O))=O"]',
    'SN1_SMILES': 'O/C=C\\CCCCCCCCCCCCCC', 'LPP_ORIGIN': 'PE(P-16:0/18:4)', 'SN1_ABBR': 'P-16:0',
    'SN_JSON': '{"SN1": "UNMOD", "SN2": "OCP"}', 'LPP_CLASS': 'PE', 'SN1_FRAGS': '[""]',
    'SN1_JSON': '{"C": 16, "KETO": 0, "OH": 0, "OAP": 0, "OCP": 0, "COOH": 0, "DB": 0, "LINK_TYPE": "P-", "CHO": 0}',
    'LM_ID': 'PE(P-16:0/15:3[3xDB,1xOH]<COOH@C15>)',
    'SN2_JSON': '{"C": 15, "KETO": 0, "OH": 1, "OAP": 3, "OCP": 1, "COOH": 1, "DB": 3, "LINK_TYPE": "", "CHO": 0}',
    'LPP_SMILES': 'OP(OCCN)(OCC([H])(OC(CCCC/C=C/C/C=C/C(O)/C=C/CC(O)=O)=O)CO/C=C\\CCCCCCCCCCCCCC)=O',
    'SN2_SMILES': 'OC(CCCC/C=C/C/C=C/C(O)/C=C/CC(O)=O)=O', 'SN2_ABBR': '15:3[3xDB,1xOH]<COOH@C15>'}

lpp_dct2 = {
    'LPP_FRAG': '["OP(OCCN)(OCC([H])(OC(CC)=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O))=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/)=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O))=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C)=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C/C=C/CC(O)=O)=O)COC(CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CC)=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O))=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/)=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O))=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C)=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C/C=C/CC(O)=O)=O)COC(CCCCCCC(O))=O)=O", "OP(OCCN)(OCC([H])(OC(CC)=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O))=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/)=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O))=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C)=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O", "OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C/C=C/CC(O)=O)=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O"]',
    'SN2_FRAGS': '["OC(CC)=O", "OC(CCC(O))=O", "OC(CCC(O)/C=C/)=O", "OC(CCC(O)/C=C/C(O))=O", "OC(CCC(O)/C=C/C(O)/C=C/C)=O"]',
    'SN1_SMILES': 'OC(CCCCCCC(O)/C=C/CCCCCC)=O', 'LPP_ORIGIN': 'PE(16:1/20:4)', 'SN1_ABBR': '16:1[1xDB,1xOH]',
    'SN_JSON': '{"SN1": "OAP", "SN2": "OCP"}', 'LPP_CLASS': 'PE', 'SN1_FRAGS': '["OC(CCCCCC)=O","OC(CCCCCCC(O))=O"]',
    'SN1_JSON': '{"C": 16, "KETO": 0, "OH": 1, "OAP": 1, "OCP": 0, "COOH": 0, "DB": 1, "LINK_TYPE": "", "CHO": 0}',
    'LM_ID': 'PE(16:1[1xDB,1xOH]/14:3[3xDB,2xOH]<COOH@C14>)',
    'SN2_JSON': '{"C": 14, "KETO": 0, "OH": 2, "OAP": 3, "OCP": 1, "COOH": 1, "DB": 3, "LINK_TYPE": "", "CHO": 0}',
    'LPP_SMILES': 'OP(OCCN)(OCC([H])(OC(CCC(O)/C=C/C(O)/C=C/C/C=C/CC(O)=O)=O)COC(CCCCCCC(O)/C=C/CCCCCC)=O)=O',
    'SN2_SMILES': 'OC(CCC(O)/C=C/C(O)/C=C/C/C=C/CC(O)=O)=O', 'SN2_ABBR': '14:3[3xDB,2xOH]<COOH@C14>'}

lpp_dct3 = {'LPP_FRAG': '["OP(OCCN)(OCC([H])(OC(CCCCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCC/C=C\C(O)CCCCCC)=O)=O"]',
            'SN2_FRAGS': '[""]', 'SN1_SMILES': 'OC(CCCCCCCCCCCCCCCCCCC)=O', 'LPP_ORIGIN': 'PE(16:0/20:0)', 'SN1_ABBR': '0:0',
            'SN_JSON': '{"SN1": "LYSO", "SN2": "UNMOD"}', 'LPP_CLASS': 'PE', 'SN1_FRAGS': '[""]',
            'SN1_JSON': '{"C": 0, "KETO": 0, "OH": 0, "OAP": 0, "OCP": 0, "COOH": 0, "DB": 0, "LINK_TYPE": "", "CHO": 0}',
            'LM_ID': 'PE(0:0/18:1)',
            'SN2_JSON': '{"C": 18, "KETO": 0, "OH": 1, "OAP": 1, "OCP": 0, "COOH": 0, "DB": 1, "LINK_TYPE": "", "CHO": 0}',
            'LPP_SMILES': 'OP(OCCN)(OCC([H])(OC(CCCCCCCCCCCCCCCCCCC)=O)COC(CCCCCCCC/C=C\C(O)CCCCCC)=O)=O',
            'SN2_SMILES': 'OC(CCCCCCCC/C=C\C(O)CCCCCC)=O', 'SN2_ABBR': '20:0'}

lpp_dct4 = {'LPP_FRAG': '["OP(OCCN)(OCC([H])(OC(CCCCCCCCCCCCCCCCCCC)=O)COC(CC/C=C\C/C=C\C/C=C\C(O)CCCCCC)=O)=O"]',
            'SN2_FRAGS': '[""]', 'SN1_SMILES': 'O', 'LPP_ORIGIN': 'PE(16:0/20:0)', 'SN1_ABBR': '0:0',
            'SN_JSON': '{"SN1": "LYSO", "SN2": "UNMOD"}', 'LPP_CLASS': 'PE', 'SN1_FRAGS': '[""]',
            'SN1_JSON': '{"C": 0, "KETO": 0, "OH": 0, "OAP": 0, "OCP": 0, "COOH": 0, "DB": 0, "LINK_TYPE": "LYSO", "CHO": 0}',
            'LM_ID': 'PE(0:0/18:3)',
            'SN2_JSON': '{"C": 18, "KETO": 0, "OH": 1, "OAP": 1, "OCP": 0, "COOH": 0, "DB": 3, "LINK_TYPE": "", "CHO": 0}',
            'LPP_SMILES': 'OP(OCCN)(OCC([H])(OC(CCCCCCCCCCCCCCCCCCC)=O)COC(CC/C=C\C/C=C\C/C=C\C(O)CCCCCC)=O)=O',
            'SN2_SMILES': 'OC(CC/C=C\C/C=C\C/C=C\C(O)CCCCCC)=O', 'SN2_ABBR': '20:0'}

frag_score = r'D:\theolpp\TheoFragPatterns_csv\ion_scores_df.xlsx'

_lpp_lst = [lpp_dct1, lpp_dct2, lpp_dct3, lpp_dct4]
a = SNFrag('PE', frag_score)
for l in _lpp_lst:
    a.calc_frags(l)
