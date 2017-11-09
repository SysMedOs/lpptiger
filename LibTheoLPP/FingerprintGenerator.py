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

import json
import re

import pandas as pd


class FingerprintGen(object):
    def __init__(self, pl_hg_cfg):
        self.pl_hg_cfg_df = pd.read_excel(pl_hg_cfg)

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
        # PA H3O4P, PC C5H14NO4P
        self.pl_hg_mz_dct = {'PA': 3 * 1.0078250321 + 4 * 15.9949146221 + 30.97376151,
                             'PC': 5 * 12 + 14 * 1.0078250321 + 14.0030740052 + 4 * 15.9949146221 + 30.97376151,
                             'PE': 2 * 12 + 8 * 1.0078250321 + 14.0030740052 + 4 * 15.9949146221 + 30.97376151,
                             'PG': 3 * 12 + 9 * 1.0078250321 + 6 * 15.9949146221 + 30.97376151,
                             'PS': 3 * 12 + 8 * 1.0078250321 + 14.0030740052 + 6 * 15.9949146221 + 30.97376151}

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

        if charge == '' or charge == 'Neutral':
            _ion_elem_dct = ion_dct

        else:
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
        return ion_mz

    def get_fingerprint(self, lpp_info_dct):

        pl_class = lpp_info_dct['LPP_CLASS']
        if pl_class in self.pl_hg_mz_dct.keys():
            if pl_class == 'PS':
                pl_hg_mz = self.pl_hg_mz_dct['PA']
            elif pl_class == 'PC':
                pl_hg_mz = self.pl_hg_mz_dct['PC'] - 12.0 - 2 * 1.0078250321
            else:
                pl_hg_mz = self.pl_hg_mz_dct[pl_class]
        else:
            pl_hg_mz = self.pl_hg_mz_dct['PA']

        sn1_dct = json.loads(lpp_info_dct['SN1_JSON'])
        sn2_dct = json.loads(lpp_info_dct['SN2_JSON'])

        sn1_formula = lpp_info_dct['SN1_FORMULA']
        sn2_formula = lpp_info_dct['SN2_FORMULA']
        sn1_formula_dct = self.parse_formula(sn1_formula)
        sn2_formula_dct = self.parse_formula(sn2_formula)
        sn1_n_mz = self.formula_to_mz(sn1_formula_dct, charge='')
        sn2_n_mz = self.formula_to_mz(sn2_formula_dct, charge='')
        sn1_c_mz = self.formula_to_mz(sn1_formula_dct, charge='[M-H]-')
        sn2_c_mz = self.formula_to_mz(sn2_formula_dct, charge='[M-H]-')

        nl_water_mz = 2 * 1.0078250321 + 15.9949146221
        nl_co2_mz = 12.0 + 2 * 15.9949146221

        pl_hg_df = self.pl_hg_cfg_df[self.pl_hg_cfg_df['CLASS'] == pl_class]
        pl_frag_df = pl_hg_df[pl_hg_df['TYPE'] == 'FRAG']
        pl_frag_mz_lst = pl_frag_df['EXACTMASS'].tolist()

        sn1_frag_mz_lst = [sn1_c_mz]
        sn2_frag_mz_lst = [sn2_c_mz]
        sn3_frag_mz_lst = pl_frag_mz_lst

        glycerol_mz = 3 * 12.0 + 2 * 1.0078250321
        sn1_nl_mz_lst = [0, nl_water_mz, sn1_n_mz]
        sn2_nl_mz_lst = [0, nl_water_mz, sn2_n_mz]
        sn3_nl_mz_lst = [pl_hg_mz]

        if 'OH' in sn1_dct.keys():
            if sn1_dct['OH'] == 1:
                sn1_frag_mz_lst.append(sn1_c_mz - nl_water_mz)
                sn1_nl_mz_lst.append(sn1_n_mz - nl_water_mz)
            if sn1_dct['OH'] > 1:
                sn1_mod_count_lst = range(1, sn1_dct['OH'] + 1)
                for mod_c in sn1_mod_count_lst:
                    sn1_frag_mz_lst.append(sn1_c_mz - mod_c * nl_water_mz)
                    sn1_nl_mz_lst.append(sn1_n_mz - mod_c * nl_water_mz)

        if 'OH' in sn2_dct.keys():
            if sn2_dct['OH'] == 1:
                sn2_frag_mz_lst.append(sn2_c_mz - nl_water_mz)
                sn2_nl_mz_lst.append(sn2_n_mz - nl_water_mz)
            if sn2_dct['OH'] > 1:
                sn2_mod_count_lst = range(1, sn2_dct['OH'] + 1)
                for mod_c in sn2_mod_count_lst:
                    sn2_frag_mz_lst.append(sn2_c_mz - mod_c * nl_water_mz)
                    sn2_nl_mz_lst.append(sn2_n_mz - mod_c * nl_water_mz)

        if 'COOH' in sn1_dct.keys():
            if sn1_dct['COOH'] == 1:
                sn1_frag_mz_lst.append(sn1_c_mz - nl_co2_mz)
                sn1_nl_mz_lst.append(sn1_n_mz - nl_co2_mz)

        if 'COOH' in sn2_dct.keys():
            if sn2_dct['COOH'] == 1:
                sn2_frag_mz_lst.append(sn2_c_mz - nl_co2_mz)
                sn2_nl_mz_lst.append(sn2_n_mz - nl_co2_mz)

        nl_ion_mz_lst = []
        for sn3_mz in sn3_nl_mz_lst:
            for sn2_mz in sn2_nl_mz_lst:
                for sn1_mz in sn1_nl_mz_lst:
                    nl_ion_mz_lst.append(glycerol_mz + sn3_mz + sn2_mz + sn1_mz - 1.0078250321)

        fp_mz_lst = sn1_frag_mz_lst + sn2_frag_mz_lst + sn3_frag_mz_lst + nl_ion_mz_lst

        if 'MSP_JSON' in lpp_info_dct:
            msp_info_df = pd.read_json(lpp_info_dct['MSP_JSON'], orient='index')
            msp_mz_lst = msp_info_df['mz'].tolist()
            for _msp in msp_mz_lst:
                if _msp in fp_mz_lst:
                    pass
                else:
                    fp_mz_lst.append(_msp)
        else:
            pass

        fp_mz_lst = sorted(list(set(round(mz, 3) for mz in fp_mz_lst if mz >= 90)))

        # print('sn1_frag_mz_lst', sn1_frag_mz_lst)
        # print('sn2_frag_mz_lst', sn2_frag_mz_lst)
        # print('sn3_frag_mz_lst', sn3_frag_mz_lst)
        # print('nl_ion_mz_lst', sorted(nl_ion_mz_lst))
        # print('fp_mz_lst', fp_mz_lst)

        return fp_mz_lst
