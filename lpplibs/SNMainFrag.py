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

from MergeBackLPP import pl_lpp


class SNMainFrag(object):
    def __init__(self, pl_class, frag_score_list):
        self.pl_type = pl_class
        print(self.pl_type)
        try:
            pl_frag_df = pd.read_excel(frag_score_list, sheetname=pl_class)

            pl_frag_df = pl_frag_df[pl_frag_df['PL_CLASS'] == self.pl_type]
            pl_frag_df['sn2_FA'] = pl_frag_df['sn2_FA'].astype('int')
            pl_frag_df['sn2_DB'] = pl_frag_df['sn2_DB'].astype('int')
            pl_frag_df['OH'] = pl_frag_df['OH'].astype('int')
            pl_frag_df['KETO'] = pl_frag_df['KETO'].astype('int')
            self.pl_frag_df = pl_frag_df
            print('self.pl_frag_df', self.pl_frag_df.shape)
        except:
            pass

        self.charge_elem_dct = {'[M+H]+': {'H': 1}, '[M+Na]+': {'Na': 1},
                                '[M+K]+': {'Na': 1}, '[M+NH4]+': {'H': 4, 'N': 1},
                                '[M-H]-': {'H': -1}, '[M+FA-H]-': {'H': 1, 'C': 1, 'O': 2},
                                '[M+HCOO]-': {'H': 1, 'C': 1, 'O': 2},
                                '+': {'H': 1}, '-': {'H': -1}}

        self.charge_mz_dct = {'[M+H]+': 1.007825, '[M+Na]+': 22.989770, '[M+K]+': 38.963708, '[M+NH4]+': 18.034374,
                              '[M-H]-': -1.007825, '[M+FA-H]-': 44.997654, '[M+HCOO]-': 44.997654}

        # iupac '97
        # elem, num, [exact mass isotopes], [abundances each isotope]
        self.periodic_table_dct = {'H': [1, [1.0078250321, 2.0141017780], [0.999885, 0.0001157]],
                                   'D': [1, [2.0141017780], [0.0001157]],
                                   'C': [6, [12.0, 13.0033548378], [0.9893, 0.0107]],
                                   'N': [7, [14.0030740052, 15.0001088984], [0.99632, 0.00368]],
                                   'O': [8, [15.9949146221, 16.99913150, 17.9991604], [0.99757, 0.00038, 0.00205]],
                                   'Na': [11, [22.98976967], [1.0]],
                                   'P': [15, [30.97376151], [1.0]],
                                   'S': [16, [31.97207069, 32.97145850, 33.96786683, 35.96708088],
                                         [0.9493, 0.0076, 0.0429, 0.0002]],
                                   'Cl': [17, [34.96885271, 36.96590260], [0.7578, 0.2422]],
                                   'K': [19, [38.9637069, 39.96399867, 40.96182597], [0.932581, 0.000117, 0.067302]],
                                   'Ca': [20, [39.9625912, 41.9586183, 42.9587668, 43.9554811, 45.9536928, 47.952534],
                                          [0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187]],
                                   'Fe': [26, [53.9396148, 55.9349421, 56.9353987, 57.9332805],
                                          [0.05845, 0.91754, 0.02119, 0.00282]],
                                   'Cu': [29, [62.9296011, 64.9277937], [0.6917, 0.3083]]}

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
            print('====>>>>>sn1 and sn2 both modified!! currently not supported!! ====>>>>> skip now !!')

        elif _sn1_mod_dct['LINK_TYPE'] in ['O-', 'P-'] or _sn2_mod_dct['LINK_TYPE'] in ['O-', 'P-']:
            print('====>>>>>sn1 or sn2 O-/P- link!! currently not supported!! ====>>>>> skip now !!')
        else:

            # get fragment scores in pandas series
            _score_se = self.get_score_type(_sn2_mod_dct, _lpp_type)
            # print(_score_se)

            if _sn1_mod_dct['LINK_TYPE'] == 'LYSO' or _sn2_mod_dct['LINK_TYPE'] == 'LYSO':
                print('>>>Skip lyso species now >>>>>>')
                pass
            else:
                # if _lpp_type == 'PC':
                #     _frag_dct['[M+FA-H]-'] = self.calc_mz(_lpp_full_smi, score=int(_score_se.loc['[M+FA-H]-']),
                #                                           charge='[M+FA-H]-')
                #     # _frag_dct['[M-H]-'] = self.calc_mz(_lpp_full_smi, score=int(_score_se.loc['[M-H]-']))
                # else:
                _pr_info = self.calc_mz(_lpp_full_smi, score=int(_score_se.loc['[M-H]-']))
                if _pr_info[1] > 0:
                    _frag_dct['[M-H]-'] = _pr_info
                else:
                    pass

                if '[sn1-H]-' in frag_type_lst:

                    _frag_dct['[sn1-H]-'] = self.calc_mz(_sn1_smi, score=int(_score_se.loc['[sn1-H]-']))

                if '[sn2-H]-' in frag_type_lst:
                    _frag_dct['[sn2-H]-'] = self.calc_mz(_sn2_smi, score=int(_score_se.loc['[sn2-H]-']))

                # for the rest part of modifications
                # change score of sn1 and sn2 if sn1 is mod but sn2 unmod
                if int(_sn1_mod_dct['OAP']) + int(_sn1_mod_dct['OCP']) > 0 and \
                                        int(_sn2_mod_dct['OAP']) + int(_sn2_mod_dct['OCP']) == 0:
                    _score_se = _score_se.rename({'[sn1-H]-': '[sn2-H]-', '[sn2-H]-': '[sn1-H]-'})
                    if '[sn1-H2O-H]-' in frag_type_lst and '[sn2-H2O-H]-':
                        _score_se = _score_se.rename({'[sn1-H2O-H]-': '[sn2-H2O-H]-', '[sn2-H2O-H]-': '[sn1-H2O-H]-'})
                    if '[sn1-CO2-H]-' in frag_type_lst and '[sn2-CO2-H]-':
                        _score_se = _score_se.rename({'[sn1-CO2-H]-': '[sn2-CO2-H]-', '[sn2-CO2-H]-': '[sn1-CO2-H]-'})

                    _sn_switch_dct ={}
                    for _ion_typ in frag_type_lst:
                        if _ion_typ in ['[sn1-H]-', '[sn2-H]-',
                                        '[sn1-H2O-H]-', '[sn2-H2O-H]-',
                                        '[sn1-CO2-H]-', '[sn2-CO2-H]-']:
                            pass
                        else:
                            sn1_chk = re.compile(r'(\[.*)([sS][nN]1)(.*[-])')
                            sn2_chk = re.compile(r'(\[.*)([sS][nN]2)(.*[-])')
                            if sn1_chk.match(_ion_typ):
                                _sn_typ = sn1_chk.sub(r'\1sn2\3', _ion_typ)
                                _sn_switch_dct[_ion_typ] = _sn_typ
                            elif sn2_chk.match(_ion_typ):
                                _sn_typ = sn2_chk.sub(r'\1sn1\3', _ion_typ)
                                _sn_switch_dct[_ion_typ] = _sn_typ

                    if len(_sn_switch_dct.keys()) > 0:
                        _score_se = _score_se.rename(_sn_switch_dct)

                    pre_frag_type_lst = _score_se.index.tolist()
                    frag_type_lst = []
                    for _typ in pre_frag_type_lst:
                        if _typ[-1] == '-':
                            frag_type_lst.append(_typ)
                    print('>>> sn1 mod & sn2 unmod >>>')

                for _ion_typ in frag_type_lst:
                    if _ion_typ in ['[M-H]-', '[M+FA-H]-', '[sn1-H]-', '[sn2-H]-']:
                        pass
                    else:
                        # print(_score_se.loc[_ion_typ])

                        if int(_score_se.loc[_ion_typ]) > 0:
                            _frag_dct[_ion_typ] = self.calc_mz(elem_info=None, mod=_ion_typ,
                                                               score=int(_score_se.loc[_ion_typ]),
                                                               charge='[M-H]-', lpp_info_dct=lpp_info_dct)
                        else:
                            pass

        msp_json = json.dumps(_frag_dct)
        return msp_json

    def get_score_type(self, sn_mod_dct, lpp_class):
        """
        This is to decide the possible fragments and corresponding intensities.
        The type and number of each modification will be checked.
        The corresponding row in the ion_scores_df will be selected and given to the fragments.
        :param sn_mod_dct:
        :param lpp_class:
        :return:
        """

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

        if sn_oap == 0 and sn_ocp == 0:
            _lpp_r = self.pl_frag_df[self.pl_frag_df['TYPE'] == 'UNMOD']
            return _lpp_r.iloc[0, :]

        if sn_oap >= 1:
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

    def calc_mz(self, elem_info, mod=None, score=0, charge='[M-H]-', lpp_info_dct=None):

        if charge in self.charge_mz_dct.keys() and charge in self.charge_elem_dct.keys():
            pass
        else:
            charge = '[M-H]-'

        if isinstance(elem_info, str):
            # test if elem_info is smiles code
            try:
                _mol = Chem.MolFromSmiles(elem_info)
                AllChem.Compute2DCoords(_mol)
                # _exactmass = rdMolDescriptors.CalcExactMolWt(_mol)
                _formula = rdMolDescriptors.CalcMolFormula(_mol)
                _elem_dct = self.parse_formula(_formula)
            except:
                _elem_dct = self.parse_formula(elem_info)
        elif isinstance(elem_info, dict):
            _elem_dct = elem_info.copy()
        else:
            _elem_dct = {}

        if mod is not None or mod != '':
            if _elem_dct is None or _elem_dct == {}:
                _elem_dct = self.get_mod_elem(elem_dct=None, mod=mod, lpp_info_dct=lpp_info_dct)
            else:
                _elem_dct = self.get_mod_elem(elem_dct=_elem_dct, mod=mod, lpp_info_dct=lpp_info_dct)

        ion_mz, _ion_elem_dct = self.formula_to_mz(_elem_dct, charge=charge)

        elem_order_lst = ['C', 'H', 'N', 'O', 'P', 'S', 'Na', 'K']
        _ion_elem = ''
        for _e in elem_order_lst:
            if _e in _ion_elem_dct.keys():
                _ion_elem += _e
                if _ion_elem_dct[_e] > 1:
                    _ion_elem += str(_ion_elem_dct[_e])
                else:
                    pass
        if charge in ['[M+H]+', '[M+Na]+', '[M+K]+', '[M+NH4]+']:
            _ion_elem += '+'
        elif charge in ['[M-H]-', '[M+FA-H]-', '[M+HCOO]-']:
            _ion_elem += '-'

        # charged_info = '|'.join([frag_type, _ion_elem])

        ion_info = (round(ion_mz, 4), score, _ion_elem)

        return ion_info

    def calc_mod_mz(self, origin_info, mod=None, score=0, charge='[M-H]-'):

        _origin_elem = origin_info[2]
        # _origin_mz = origin_info[0]

        if charge in self.charge_elem_dct.keys():
            pass
        else:
            charge = '[M-H]-'

        mod_dct = {'-H2O': {'H': -2, 'O': -1}, '+H2O': {'H': 2, 'O': 1}, '-CO2': {'C': -1, 'O': -2}}

        if re.match(r'.*[-]H2O.*', mod):
            _mod_elem_dct = mod_dct['-H2O']
        elif re.match(r'.*[+]H2O.*', mod):
            _mod_elem_dct = mod_dct['+H2O']
        elif re.match(r'.*[-]CO2.*', mod):
            _mod_elem_dct = mod_dct['-CO2']
        else:
            _mod_elem_dct = {}

        _origin_elem_dct = self.parse_formula(_origin_elem)
        _mod_fin_dct = _origin_elem_dct.copy()

        for _elem in _mod_elem_dct.keys():
            _mod_fin_dct[_elem] += _mod_elem_dct[_elem]

        mod_mz = self.formula_to_mz(_mod_fin_dct, charge=charge)

    @staticmethod
    def parse_formula(formula=None):

        _elem_dct = {}

        _elem_rgx = re.compile(r'[A-Z][a-z]?[0-9]{0,3}')
        _elem_bulk_lst = re.findall(_elem_rgx, formula)
        _elem_num_rgx = re.compile(r'([A-Z][a-z]?)([0-9]{0,3})')

        for _elem in _elem_bulk_lst:
            _elem_checker = re.match(_elem_num_rgx, _elem)
            _elem_lst = _elem_checker.groups()
            try:
                _elem_dct[_elem_lst[0]] = int(_elem_lst[1])
            # e.g. only 1 N or P will get '' for _elem_lst[1]
            except ValueError:
                _elem_dct[_elem_lst[0]] = 1

        return _elem_dct

    def formula_to_mz(self, elem_info, charge='[M-H]-'):

        if isinstance(elem_info, str):
            ion_dct = self.parse_formula(elem_info)
        elif isinstance(elem_info, dict):
            ion_dct = elem_info.copy()
        else:
            ion_dct = {}

        if charge in self.charge_elem_dct.keys():
            charge_dct = self.charge_elem_dct[charge]
        else:
            # set as '[M-H]-'
            charge_dct = {'H': -1}

        # get sum keys form both dict
        _ion_elem_dct = {}
        _charged_keys_lst = set(sum([ion_dct.keys(), charge_dct.keys()], []))
        for _key in _charged_keys_lst:
            _ion_elem_dct[_key] = ion_dct.get(_key, 0) + charge_dct.get(_key, 0)

        ion_mz = 0.0

        for _elem in _ion_elem_dct.keys():
            if _elem in self.periodic_table_dct.keys():
                ion_mz += _ion_elem_dct[_elem] * self.periodic_table_dct[_elem][1][0]
            else:
                print('!!Unsupported elements!!')
                break

        ion_mz = round(ion_mz, 4)
        return ion_mz, _ion_elem_dct

    def get_mod_elem(self, elem_dct=None, mod=None, lpp_info_dct=None):

        mod_dct = {'-H2O': {'H': -2, 'O': -1}, '+H2O': {'H': 2, 'O': 1},
                   '-CO2': {'C': -1, 'O': -2}, '+FA': {'H': 2, 'C': 1, 'O': 2},
                   '-C3H9N': {'C': -3, 'H': -9, 'N': -1},
                   '-C3H5NO2': {'C': -3, 'O': -2, 'H': -5, 'N': -1},
                   '-CH3COOH': {'C': -2, 'O': -2, 'H': -4},
                   '-CH2': {'C': -1, 'H': -2}}

        # get the formula as dict
        if elem_dct is None:
            if lpp_info_dct is None:
                _elem_dct = {}
            else:
                _lpp_type = lpp_info_dct['LPP_CLASS']
                _lpp_full_smi = lpp_info_dct['LPP_SMILES']
                _sn1_smi = lpp_info_dct['SN1_SMILES']
                _sn2_smi = lpp_info_dct['SN2_SMILES']
                _lyso_smi = 'O'
                if re.match(r'\[M-[sS][nN][1].*', mod):
                    _frag_smi = pl_lpp(_lpp_type, sn1=_lyso_smi, sn2=_sn2_smi)
                elif re.match(r'\[M-[sS][nN][2].*', mod):
                    _frag_smi = pl_lpp(_lpp_type, sn1=_sn1_smi, sn2=_lyso_smi)
                elif re.match(r'\[[sS][nN][1].*', mod):
                    _frag_smi = _sn1_smi
                elif re.match(r'\[[sS][nN][2].*', mod):
                    _frag_smi = _sn2_smi
                elif re.match(r'\[M[+-][HCFA].*', mod):
                    _frag_smi = _lpp_full_smi
                else:
                    _frag_smi = ''

                _mol = Chem.MolFromSmiles(_frag_smi)
                AllChem.Compute2DCoords(_mol)
                _formula = rdMolDescriptors.CalcMolFormula(_mol)
                _elem_dct = self.parse_formula(_formula)
        else:
            _elem_dct = elem_dct.copy()

        # print('mod', mod)
        _mod_sum_elem_dct = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
        if mod is None or mod == '':
            _mod_sum_elem_dct = {}
        else:
            chk0 = re.compile(r'.*-[sS][nN][12].*')
            chk1 = re.compile(r'.*[-]H2O.*')
            chk2 = re.compile(r'.*[+]H2O.*')
            chk3 = re.compile(r'.*[-]CO2.*')
            chk4 = re.compile(r'.*[-]C3H9N.*')
            chk5 = re.compile(r'.*[-]CH3COOH.*')
            chk6 = re.compile(r'.*[-]CH2.*')
            chk7 = re.compile(r'.*[-]C3H5NO2.*')

            chk9 = re.compile(r'.*[+]FA.*')
            chk10 = re.compile(r'(P[ACEGSI]4?P?_)(.*)([+-])')

            if chk0.match(mod):
                _mod_elem_dct = mod_dct['-H2O']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]

            if chk1.match(mod):
                _mod_elem_dct = mod_dct['-H2O']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            # for [M-sn+H2O-H]-, the H2O was add already
            # elif chk2.match(mod):
            #     _mod_elem_dct = mod_dct['+H2O']
            if chk2.match(mod):
                _mod_elem_dct = mod_dct['+H2O']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk3.match(mod):
                _mod_elem_dct = mod_dct['-CO2']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk4.match(mod):
                _mod_elem_dct = mod_dct['-C3H9N']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk5.match(mod):
                _mod_elem_dct = mod_dct['-CH3COOH']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk6.match(mod):
                _mod_elem_dct = mod_dct['-CH2']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk7.match(mod):
                _mod_elem_dct = mod_dct['-C3H5NO2']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]

            if chk9.match(mod):
                _mod_elem_dct = mod_dct['+FA']
                for _key in _mod_elem_dct.keys():
                    _mod_sum_elem_dct[_key] += _mod_elem_dct[_key]
            if chk10.match(mod):
                m7 = chk10.match(mod)
                m7_lst = m7.groups()
                m7_elem = m7_lst[1]
                _charge = m7_lst[2]
                _mod_sum_elem_dct = self.parse_formula(formula=m7_elem)
                # naturalize the formula
                if _charge in self.charge_elem_dct.keys():
                    _mod_sum_elem_dct['H'] += -1 * self.charge_elem_dct[_charge]['H']
                # print('found PE HG ------------------------------------------------------>>>>>>>>>>>>>>>')
            else:
                pass

        # get sum keys form both dict
        _frag_elem_dct = {}
        _charged_keys_lst = set(sum([_elem_dct.keys(), _mod_sum_elem_dct.keys()], []))
        for _key in _charged_keys_lst:
            _frag_elem_dct[_key] = _elem_dct.get(_key, 0) + _mod_sum_elem_dct.get(_key, 0)

        return _frag_elem_dct
