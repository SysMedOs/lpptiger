# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPsmi.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from __future__ import division

import re
# import json
import math

import numpy as np
import pandas as pd
from scipy import spatial


class ScoreGenerator:
    def __init__(self, fa_def_df, weight_df, key_frag_df, lipid_type, ion_charge='[M-H]-', ms2_precision=500e-6):
        self.fa_def_df = fa_def_df
        self.weight_dct = weight_df.to_dict()
        self.weight_type_lst = self.weight_dct.keys()
        self.target_frag_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "FRAG" and PR_CHARGE == "%s"'
                                                % (lipid_type, ion_charge))
        self.target_nl_df = key_frag_df.query(r'CLASS == "%s" and TYPE == "NL" and PR_CHARGE == "%s"'
                                              % (lipid_type, ion_charge))

        if ion_charge in ['[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+FA-H]-', '[M+OAc]-']:
            charge_mode = 'NEG'
        elif ion_charge in ['[M+H]+', '[M+NH4]+']:
            charge_mode = 'POS'
        else:
            charge_mode = 'NEG'

        self.other_frag_df = key_frag_df.query('CLASS != "%s" and TYPE == "FRAG" and CHARGE_MODE == "%s"'
                                               % (lipid_type, charge_mode))
        self.other_nl_df = key_frag_df.query('CLASS != "%s" and TYPE == "NL" and CHARGE_MODE == "%s"'
                                             % (lipid_type, charge_mode))
        self.lipid_type = lipid_type

        self.fa_uniquemz_df = pd.DataFrame(self.fa_def_df, columns=['elem', 'mass', '[M-H]-', 'NL-H2O'])
        # Keep same elem only once to speed up
        self.fa_uniquemz_df.drop_duplicates(subset='elem', inplace=True)
        self.fa_uniquemz_df.rename(columns={'mass': 'NL'}, inplace=True)
        self.fa_uniquemz_df.loc[:, 'sn'] = self.fa_uniquemz_df.loc[:, '[M-H]-']
        self.fa_uniquemz_df.loc[:, '[M-H]-_s'] = self.fa_uniquemz_df['[M-H]-'] * (1 - ms2_precision)
        self.fa_uniquemz_df.loc[:, '[M-H]-_b'] = self.fa_uniquemz_df['[M-H]-'] * (1 + ms2_precision)
        self.fa_uniquemz_df.loc[:, 'NL_s'] = self.fa_uniquemz_df['NL'] * (1 - ms2_precision)
        self.fa_uniquemz_df.loc[:, 'NL_b'] = self.fa_uniquemz_df['NL'] * (1 + ms2_precision)
        self.fa_uniquemz_df.loc[:, 'NL-H2O_s'] = self.fa_uniquemz_df['NL-H2O'] * (1 - ms2_precision)
        self.fa_uniquemz_df.loc[:, 'NL-H2O_b'] = self.fa_uniquemz_df['NL-H2O'] * (1 + ms2_precision)
        for _i_u_mz, _u_mz_r in self.fa_uniquemz_df.iterrows():
            self.fa_uniquemz_df.set_value(_i_u_mz, '[M-H]-_query', '%f <= mz <= %f' % (_u_mz_r['[M-H]-_s'],
                                                                                       _u_mz_r['[M-H]-_b']))

    @staticmethod
    def get_pr_mz(charge_type, mz_lib):

        pr_mz = 0.0

        if charge_type in ['[M-H]-', '[M+HCOO]-', '[M+FA-H]-', '[M+CH3COO]-', '[M+OAc-H]-', '[M+AcOH-H]-']:
            charge_mode = 'NEG'
            if charge_type == '[M-H]-':
                pr_mz = mz_lib
            elif charge_type in ['[M+HCOO]-', '[M+FA-H]-']:
                pr_mz = mz_lib - 46.005480  # - HCOOH
            elif charge_type == '[M+CH3COO]-':
                pr_mz = mz_lib - 60.021130  # - CH3COOH

        elif charge_type in ['[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+']:
            charge_mode = 'POS'
            if charge_type == '[M+H]+':
                pr_mz = mz_lib
            elif charge_type == '[M+Na]+':
                pr_mz = mz_lib - 22.989770 + 1.007825  # - Na + H
            elif charge_type == '[M+NH4]+':
                pr_mz = mz_lib - 17.026549  # - NH3
            elif charge_type == '[M+K]+':
                pr_mz = mz_lib - 38.963708 + 1.007825  # - K + H
        else:
            charge_mode = 'NEG'
            pr_mz = mz_lib

        return pr_mz, charge_mode

    @staticmethod
    def decode_abbr(abbr):

        pl_checker = re.compile(r'(P[ACEGSI])([(])(.*)([)])')
        pip_checker = re.compile(r'(PIP)([(])(.*)([)])')
        tg_checker = re.compile(r'(TG)([(])(.*)([)])')
        fa_checker = re.compile(r'(.*)(\d{1,2})([:])(\d{1,2})(.*)([/_])(.*)(\d{1,2})([:])(\d{1,2})(.*)')
        fa_o_checker = re.compile(r'(O-)(\d{1,2})([:])(\d)(.*)([/_])(\d{1,2})([:])(\d{1,2})(.*)')
        fa_p_checker = re.compile(r'(P-)(\d{1,2})([:])(\d)(.*)([/_])(\d{1,2})([:])(\d{1,2})(.*)')

        # Check PL Type
        _pl_typ = ''
        info_fa_typ = ''
        bulk_fa_linker = ''
        sn1_fa_abbr = ''
        lyso_fa_linker_dct = {'sn1': '', 'sn2': ''}

        if pl_checker.match(abbr):
            # print('PL')
            pl_re_chk = pl_checker.match(abbr)
            pl_typ_lst = pl_re_chk.groups()
            _pl_typ = pl_typ_lst[0]
            info_fa_typ = pl_typ_lst[2]
        if pip_checker.match(abbr):
            # print('PIP')
            pip_re_chk = pip_checker.match(abbr)
            pip_typ_lst = pip_re_chk.groups()
            _pl_typ = pip_typ_lst[0]
            info_fa_typ = pip_typ_lst[2]
        if tg_checker.match(abbr):
            # print('TG')
            tg_re_chk = tg_checker.match(abbr)
            tg_typ_lst = tg_re_chk.groups()
            _pl_typ = tg_typ_lst[0]
            info_fa_typ = tg_typ_lst[2]

        if fa_checker.match(info_fa_typ):

            if fa_o_checker.match(info_fa_typ):
                bulk_fa_linker = 'O-A-'
                lyso_fa_linker_dct = {'O': '', 'A': 'O-'}  # link of the other sn after NL of this sn
                fa_chk = fa_o_checker.match(info_fa_typ)
                info_fa_lst = fa_chk.groups()
                sn1_fa_abbr = ''.join(info_fa_lst[0:3])
                sn2_fa_abbr = ''.join(info_fa_lst[5:])
            elif fa_p_checker.match(info_fa_typ):
                bulk_fa_linker = 'P-A-'
                lyso_fa_linker_dct = {'P': '', 'A': 'P-'}  # link of the other sn after NL of this sn
                fa_chk = fa_p_checker.match(info_fa_typ)
                info_fa_lst = fa_chk.groups()
                sn2_fa_abbr = ''.join(info_fa_lst[5:])
            else:
                bulk_fa_linker = 'A-A-'
                lyso_fa_linker_dct = {'A': ''}
                fa_chk = fa_checker.match(info_fa_typ)
                info_fa_lst = fa_chk.groups()
                sn1_fa_abbr = ''.join(info_fa_lst[0:5])
                sn2_fa_abbr = ''.join(info_fa_lst[6:])

        else:
            sn1_fa_abbr = ''
            sn2_fa_abbr = ''

        lipid_info_dct = {'TYPE': _pl_typ, 'LINK': bulk_fa_linker,
                          'sn1_abbr': sn1_fa_abbr, 'sn2_abbr': sn2_fa_abbr,
                          'LYSO_LINK': lyso_fa_linker_dct}

        return lipid_info_dct

    def get_structure(self, abbr):

        print(abbr)

        lipid_info_dct = self.decode_abbr(abbr)
        pl_typ = lipid_info_dct['TYPE']
        sn1_abbr = lipid_info_dct['sn1_abbr']
        sn2_abbr = lipid_info_dct['sn2_abbr']

        fa_abbr_lst = self.fa_def_df['FA'].tolist()
        if pl_typ in ['PE', 'PA', 'PC', 'PI', 'PS', 'PG']:
            if sn1_abbr in fa_abbr_lst and sn2_abbr in fa_abbr_lst:

                rebuild_pl = ''.join([pl_typ, '(', sn1_abbr, '/', sn2_abbr, ')'])

                if rebuild_pl == abbr:
                    lipid_abbr_dct = {'Proposed_structures': rebuild_pl, 'Class': pl_typ,
                                      'sn1_abbr': sn1_abbr, 'sn2_abbr': sn2_abbr, }

                else:
                    lipid_abbr_dct = {}
            else:
                lipid_abbr_dct = {}
        else:
            lipid_abbr_dct = {}

        return lipid_abbr_dct

    def get_fa_signals(self, lipid_info_dct, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
                       ms2_threshold=100, ms2_infopeak_threshold=0.02):

        ident_df_dct = {}
        ident_checker = 0

        sn1_fa = lipid_info_dct['sn1_abbr']
        sn2_fa = lipid_info_dct['sn2_abbr']
        pl_class = lipid_info_dct['Class']

        sn1_mz = self.fa_def_df.loc[self.fa_def_df['FA'] == sn1_fa, '[M-H]-'].values[0]
        sn1_nl = self.fa_def_df.loc[self.fa_def_df['FA'] == sn1_fa, 'mass'].values[0]
        sn1_nl_w = self.fa_def_df.loc[self.fa_def_df['FA'] == sn1_fa, 'NL-H2O'].values[0]

        sn2_mz = self.fa_def_df.loc[self.fa_def_df['FA'] == sn2_fa, '[M-H]-'].values[0]
        sn2_nl = self.fa_def_df.loc[self.fa_def_df['FA'] == sn2_fa, 'mass'].values[0]
        sn2_nl_w = self.fa_def_df.loc[self.fa_def_df['FA'] == sn2_fa, 'NL-H2O'].values[0]

        sn_mz_dct = {'sn1': sn1_mz, 'sn2': sn2_mz}

        if pl_class in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
            lipid_type = 'PL'
        elif pl_class in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
            lipid_type = 'GL'
        else:
            lipid_type = 'PL'

        if lipid_type == 'PL':

            if pl_class == 'PC' and charge_type == '[M+HCOO]-':
                sn_mz_dct['[M-H]-sn1'] = mz_lib - sn1_nl_w - 60.021130  # - CH3COOH for PC
                sn_mz_dct['[M-H]-sn2'] = mz_lib - sn2_nl_w - 60.021130  # - CH3COOH for PC
                sn_mz_dct['[M-H]-sn1-H2O'] = mz_lib - sn1_nl - 60.021130  # - CH3COOH for PC
                sn_mz_dct['[M-H]-sn2-H2O'] = mz_lib - sn2_nl - 60.021130  # - CH3COOH for PC
            elif pl_class == 'PS' and charge_type == '[M-H]-':
                sn_mz_dct['[M-H]-sn1'] = mz_lib - sn1_nl_w - 87.032029  # - C3H5NO2 for PS
                sn_mz_dct['[M-H]-sn2'] = mz_lib - sn2_nl_w - 87.032029  # - C3H5NO2 for PS
                sn_mz_dct['[M-H]-sn1-H2O'] = mz_lib - sn1_nl - 87.032029  # - C3H5NO2 for PS
                sn_mz_dct['[M-H]-sn2-H2O'] = mz_lib - sn2_nl - 87.032029  # - C3H5NO2 for PS
            else:
                # Loss of FA-18, -OH remains on Glycerol back bone
                sn_mz_dct['[M-H]-sn1'] = mz_lib - sn1_nl_w
                sn_mz_dct['[M-H]-sn2'] = mz_lib - sn2_nl_w
                # Loss of FA as full acid, -OH remains on FA NL
                sn_mz_dct['[M-H]-sn1-H2O'] = mz_lib - sn1_nl
                sn_mz_dct['[M-H]-sn2-H2O'] = mz_lib - sn2_nl

            for _frag_type in ['sn1', '[M-H]-sn1', '[M-H]-sn1-H2O', 'sn2', '[M-H]-sn2', '[M-H]-sn2-H2O']:
                _frag_mz = sn_mz_dct[_frag_type]
                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
                _frag_mz_query_code = '%f <= mz <= %f' % (_frag_mz_low, _frag_mz_high)
                _frag_df = ms2_df.query(_frag_mz_query_code)

                if _frag_df.shape[0] > 0:
                    _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                    _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()

                    if _frag_df.shape[0] > 1:
                        _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                        _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                        _frag_df = _frag_i_df.copy()
                        if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                            pass
                        else:
                            _frag_df = _frag_i_df.append(_frag_ppm_df)
                        # convert df to dict
                        ident_df_dct[_frag_type] = _frag_df.iloc[0, :].to_dict()

                        if _frag_type == 'sn1':
                            ident_df_dct[_frag_type]['Proposed_structures'] = sn1_fa
                            ident_checker += 1
                        elif _frag_type == 'sn2':
                            ident_df_dct[_frag_type]['Proposed_structures'] = sn2_fa
                            ident_checker += 1
                        else:
                            ident_df_dct[_frag_type]['Proposed_structures'] = _frag_type

                        del _frag_df

        return ident_df_dct, ident_checker

    def get_fa_possibilities(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
                             ms2_threshold=100, ms2_infopeak_threshold=0.02):

        other_fa_df = pd.DataFrame()
        other_lyso_l_df = pd.DataFrame()
        other_lyso_h_df = pd.DataFrame()

        # use the max threshold from abs & relative intensity settings
        ms2_basepeak_i = ms2_df['i'].max()
        ms2_info_i = ms2_basepeak_i * ms2_infopeak_threshold
        ms2_threshold = max(ms2_threshold, ms2_info_i)

        calc_pr_mz = mz_lib

        if abbr[:2] in ['PA', 'PC', 'PE', 'PG', 'PI', 'PS', 'SM']:
            lipid_type = 'PL'
        elif abbr[:2] in ['TA', 'TG', 'DA', 'DG', 'MA', 'MG']:
            lipid_type = 'GL'
        else:
            lipid_type = 'PL'

        if lipid_type == 'PL':
            fa_chk_df = self.fa_uniquemz_df[['sn', 'NL', 'NL_s', 'NL_b', 'NL-H2O', 'NL-H2O_s', 'NL-H2O_b',
                                             '[M-H]-', '[M-H]-_query']]

            if abbr[:2] == 'PC' and charge_type == '[M+HCOO]-':
                fa_chk_df.loc[:, '[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 60.021130  # - CH3COOH for PC
                fa_chk_df.loc[:, '[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 60.021130  # - CH3COOH for PC
                for _i_pc, _r_pc in fa_chk_df.iterrows():
                    fa_chk_df.set_value(_i_pc, '[M-H]-sn_query', ('%f <= mz <= %f' %
                                                                  (calc_pr_mz - _r_pc['NL-H2O_b'] - 60.021130,
                                                                   calc_pr_mz - _r_pc['NL-H2O_s'] - 60.021130)))
                    fa_chk_df.set_value(_i_pc, '[M-H]-sn-H2O_query', ('%f <= mz <= %f' %
                                                                      (calc_pr_mz - _r_pc['NL_b'] - 60.021130,
                                                                       calc_pr_mz - _r_pc['NL_s'] - 60.021130)))
            elif abbr[:2] == 'PS' and charge_type == '[M-H]-':
                fa_chk_df.loc[:, '[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O'] - 87.032029  # - C3H5NO2 for PS
                fa_chk_df.loc[:, '[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL'] - 87.032029  # - C3H5NO2 for PS
                for _i_ps, _r_ps in fa_chk_df.iterrows():
                    fa_chk_df.set_value(_i_ps, '[M-H]-sn_query', ('%f <= mz <= %f' %
                                                                  (calc_pr_mz - _r_ps['NL-H2O_b'] - 87.032029,
                                                                   calc_pr_mz - _r_ps['NL-H2O_s'] - 87.032029)))
                    fa_chk_df.set_value(_i_ps, '[M-H]-sn-H2O_query', ('%f <= mz <= %f' %
                                                                      (calc_pr_mz - _r_ps['NL_b'] - 87.032029,
                                                                       calc_pr_mz - _r_ps['NL_s'] - 87.032029)))
            else:
                # Loss of FA-18, -OH remains on Glycerol back bone
                fa_chk_df.loc[:, '[M-H]-sn'] = calc_pr_mz - fa_chk_df['NL-H2O']
                # Loss of FA as full acid, -OH remains on FA NL
                fa_chk_df.loc[:, '[M-H]-sn-H2O'] = calc_pr_mz - fa_chk_df['NL']
                for _i_pr, _r_pr in fa_chk_df.iterrows():
                    fa_chk_df.set_value(_i_pr, '[M-H]-sn_query', ('%f <= mz <= %f' %
                                                                  (calc_pr_mz - _r_pr['NL-H2O_b'],
                                                                   calc_pr_mz - _r_pr['NL-H2O_s'])))
                    fa_chk_df.set_value(_i_pr, '[M-H]-sn-H2O_query', ('%f <= mz <= %f' %
                                                                      (calc_pr_mz - _r_pr['NL_b'],
                                                                       calc_pr_mz - _r_pr['NL_s'])))

            for _i, _fa_se in fa_chk_df.iterrows():

                for _frag_type in ['sn', '[M-H]-sn', '[M-H]-sn-H2O']:
                    _frag_mz = _fa_se[_frag_type]
                    if _frag_type == 'sn':
                        _frag_df = ms2_df.query(_fa_se['[M-H]-_query'])
                    elif _frag_type == '[M-H]-sn':
                        _frag_df = ms2_df.query(_fa_se['[M-H]-sn_query'])
                    elif _frag_type == '[M-H]-sn-H2O':
                        _frag_df = ms2_df.query(_fa_se['[M-H]-sn-H2O_query'])

                    if _frag_df.shape[0] > 0:
                        _frag_df.loc[:, 'ppm'] = 1e6 * (_frag_df['mz'] - _frag_mz) / _frag_mz
                        _frag_df.loc[:, 'ppm_abs'] = _frag_df['ppm'].abs()

                        if _frag_df.shape[0] > 1:
                            _frag_i_df = _frag_df.sort_values(by='i', ascending=False).head(1)
                            _frag_ppm_df = _frag_df.sort_values(by='ppm_abs').head(1)
                            _frag_df = _frag_i_df.copy()
                            if _frag_ppm_df['i'].tolist() == _frag_i_df['i'].tolist():
                                pass
                            else:
                                _frag_df = _frag_i_df.append(_frag_ppm_df)
                        if _frag_type == 'sn':
                            other_fa_df = other_fa_df.append(_frag_df)
                        elif _frag_type == '[M-H]-sn':
                            other_lyso_h_df = other_lyso_h_df.append(_frag_df)
                        elif _frag_type == '[M-H]-sn-H2O':
                            other_lyso_l_df = other_lyso_l_df.append(_frag_df)

        # format the output DataFrame
        if other_fa_df.shape[0] > 0:
            other_fa_df = other_fa_df.query('i > %f' % ms2_threshold)
            other_fa_df = other_fa_df.loc[:, ['mz', 'i', 'ppm', 'ppm_abs']]
            other_fa_df = other_fa_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            other_fa_df = other_fa_df.sort_values(by='i', ascending=False).head(10).reset_index(drop=True)

        if other_lyso_h_df.shape[0] > 0:
            other_lyso_h_df = other_lyso_h_df.query('i > %f' % ms2_threshold)
            other_lyso_h_df = other_lyso_h_df.loc[:, ['mz', 'i', 'ppm', 'ppm_abs']]
            other_lyso_h_df = other_lyso_h_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            other_lyso_h_df = other_lyso_h_df.sort_values(by='i', ascending=False).head(5).reset_index(drop=True)

        if other_lyso_l_df.shape[0] > 0:
            other_lyso_l_df = other_lyso_l_df.query('i > %f' % ms2_threshold)
            other_lyso_l_df = other_lyso_l_df.loc[:, ['mz', 'i', 'ppm', 'ppm_abs']]
            other_lyso_l_df = other_lyso_l_df.sort_values(by=['i', 'ppm_abs'], ascending=[False, True])
            other_lyso_l_df = other_lyso_l_df.sort_values(by='i', ascending=False).head(5).reset_index(drop=True)

        other_signals_dct = {'other_fa_df': other_fa_df,
                             'other_lyso_h_df': other_lyso_h_df, 'other_lyso_l_df': other_lyso_l_df}

        return other_signals_dct

    def get_rankscore(self, abbr, charge_type, mz_lib, ms2_df, ms2_precision=500e-6,
                      ms2_threshold=100, ms2_infopeak_threshold=0.02):

        lipid_info_dct = self.get_structure(abbr)

        ident_signals_dct, ident_checker = self.get_fa_signals(lipid_info_dct, charge_type, mz_lib, ms2_df,
                                                               ms2_precision=ms2_precision,
                                                               ms2_threshold=ms2_threshold,
                                                               ms2_infopeak_threshold=ms2_infopeak_threshold
                                                               )
        # need at least one FA from signal checker
        if ident_checker > 0:
            ident_signal_lst = ident_signals_dct.keys()
            if len(ident_signal_lst) > 0:

                other_signals_dct = self.get_fa_possibilities(abbr, charge_type, mz_lib, ms2_df,
                                                              ms2_precision=ms2_precision, ms2_threshold=ms2_threshold,
                                                              ms2_infopeak_threshold=ms2_infopeak_threshold)

                other_fa_df = other_signals_dct['other_fa_df']
                other_lyso_l_df = other_signals_dct['other_lyso_l_df']
                other_lyso_h_df = other_signals_dct['other_lyso_h_df']

                rank_score = 0

                matched_fa_df = pd.DataFrame()
                matched_lyso_df = pd.DataFrame()
                other_signals_df = pd.DataFrame()

                signal_mz_lst = []

                for signal in ident_signal_lst:

                    if signal in self.weight_type_lst:

                        other_signals_df = pd.DataFrame()

                        signal_dct = ident_signals_dct[signal]
                        signal_mz = signal_dct['mz']
                        signal_mz_lst.append(signal_mz)

                        if signal in ['sn1', 'sn2']:
                            try:
                                sig_idx = other_fa_df['mz'].tolist().index(signal_mz)
                            except ValueError:
                                try:
                                    sig_idx = other_fa_df['mz'].tolist().index(round(signal_mz, 6))
                                except ValueError:
                                    sig_idx = 100
                            if sig_idx <= 10:
                                rank_score += self.weight_dct[signal] * (10 - sig_idx) / 10
                                matched_fa_df = matched_fa_df.append(pd.DataFrame(data=signal_dct, index=[sig_idx]))
                                other_fa_df = other_fa_df.drop(other_fa_df.index[sig_idx])
                        if signal in ['[M-H]-sn1', '[M-H]-sn2']:
                            try:
                                sig_idx = other_lyso_h_df['mz'].tolist().index(signal_mz)
                                rank_score += self.weight_dct[signal] * (10 - sig_idx) / 10
                                matched_lyso_df = matched_lyso_df.append(pd.DataFrame(data=signal_dct, index=[sig_idx]))
                                other_lyso_h_df = other_lyso_h_df.drop(other_lyso_h_df.index[sig_idx])
                            except ValueError:
                                pass
                        if signal in ['[M-H]-sn1-H2O', '[M-H]-sn2-H2O']:
                            try:
                                sig_idx = other_lyso_l_df['mz'].tolist().index(signal_mz)
                                rank_score += self.weight_dct[signal] * (10 - sig_idx) / 10
                                matched_lyso_df = matched_lyso_df.append(pd.DataFrame(data=signal_dct, index=[sig_idx]))
                                other_lyso_l_df = other_lyso_l_df.drop(other_lyso_l_df.index[sig_idx])
                            except ValueError:
                                pass

                print('matched_fa_df')
                print(matched_fa_df)
                print('matched_lyso_df')
                print(matched_lyso_df)

                matched_checker = 1
                if other_fa_df.shape[0] > 0:
                    other_signals_df = other_signals_df.append(other_fa_df)
                if other_lyso_l_df.shape[0] > 0:
                    other_signals_df = other_signals_df.append(other_lyso_l_df)
                if other_lyso_h_df.shape[0] > 0:
                    other_signals_df = other_signals_df.append(other_lyso_h_df)

                other_signals_df = other_signals_df.reset_index(drop=True)

                other_drop_lst = []
                if other_signals_df.shape[0] > 0:
                    ident_signal_chk_lst = [round(s, 1) for s in signal_mz_lst]
                    other_signal_lst = other_signals_df['mz'].tolist()
                    other_signal_chk_lst = [round(s, 1) for s in other_signal_lst]
                    for sig_mz in ident_signal_chk_lst:
                        if sig_mz in other_signal_chk_lst:
                            other_sig_idx = other_signal_chk_lst.index(sig_mz)
                            other_drop_lst.append(other_sig_idx)

                other_signals_df = other_signals_df.drop(other_signals_df.index[other_drop_lst])

                rank_score = round(rank_score, 1)

                match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': rank_score,
                                  'OTHER_SIGNALS_INFO': other_signals_df,
                                  'MATCHED_FA_INFO': matched_fa_df, 'MATCHED_LYSO_INFO': matched_lyso_df}

            else:
                print('!! No structure related signals found !!')
                matched_checker = 0
                empty_df = pd.DataFrame()
                match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0,
                                  'OTHER_SIGNALS_INFO': empty_df,
                                  'MATCHED_FA_INFO': empty_df, 'MATCHED_LYSO_INFO': empty_df}
        else:
            print('!! No structure related signals found !!')
            matched_checker = 0
            empty_df = pd.DataFrame()
            match_info_dct = {'MATCH_INFO': matched_checker, 'Rank_score': 0.0,
                              'OTHER_SIGNALS_INFO': empty_df,
                              'MATCHED_FA_INFO': empty_df, 'MATCHED_LYSO_INFO': empty_df}

        return match_info_dct, matched_checker

    @staticmethod
    def get_cosine_score(msp_df, ms2_df, ms2_precision=500e-6, ms2_threshold=100, ms2_infopeak_threshold=0.02):

        lib_mz_lst = msp_df['mz'].tolist()

        ms2_basepeak_i = ms2_df['i'].max()
        ms2_info_i = ms2_basepeak_i * ms2_infopeak_threshold
        ms2_threshold = max(ms2_threshold, ms2_info_i)

        obs_score_df = pd.DataFrame()
        for mz in lib_mz_lst:
            mz_l = mz * (1 - ms2_precision)
            mz_h = mz * (1 + ms2_precision)

            tmp_df = ms2_df.query('%f <= mz <= %f and i >= %f ' % (mz_l, mz_h, ms2_threshold))

            # print(tmp_df)

            if tmp_df.shape[0] == 1:
                obs_score_df = obs_score_df.append(tmp_df)
            elif tmp_df.shape[0] > 1:
                tmp_df = tmp_df.sort_values(by='i', ascending=False)
                obs_score_df = obs_score_df.append(tmp_df.head(1))
            else:
                tmp_df = pd.DataFrame(data={'mz': [0.0], 'i': [0.0]})
                obs_score_df = obs_score_df.append(tmp_df)

        obs_i_max = obs_score_df['i'].max()
        # calc minus i for plot
        msp_df.loc[:, 'rev_abs_i'] = msp_df['i'] * -0.001 * obs_i_max

        obs_ar = obs_score_df.as_matrix()
        lib_ar = np.column_stack((lib_mz_lst, msp_df['i'].tolist()))

        # reduce matrix to 1D array
        obs_flat = np.hstack(obs_ar.T)
        lib_flat = np.hstack(lib_ar.T)

        cosine_score = 100 * (1 - spatial.distance.cosine(obs_flat, lib_flat))
        cosine_score = round(cosine_score, 1)
        if cosine_score > 0:
            pass
        else:
            cosine_score = 0
        print('msp_df')
        print(msp_df)
        print('obs_msp_df')
        print(obs_score_df)

        return cosine_score, msp_df, obs_score_df

    @staticmethod
    def get_fingerprint_score(fingerprint_lst, ms2_df, ms2_precision=500e-6,
                              ms2_threshold=100, ms2_infopeak_threshold=0.02):

        ms2_basepeak_i = ms2_df['i'].max()
        ms2_info_i = ms2_basepeak_i * ms2_infopeak_threshold
        ms2_threshold = max(ms2_threshold, ms2_info_i)

        obs_fp_lst = []
        missed_fp_lst = []

        fp_lib_lst = []
        fp_obs_lst = []

        obs_score_df = pd.DataFrame()
        for mz in fingerprint_lst:
            fp_lib_lst.append(1)
            mz_l = mz * (1 - ms2_precision)
            mz_h = mz * (1 + ms2_precision)

            tmp_df = ms2_df.query('%f <= mz <= %f and i >= %f ' % (mz_l, mz_h, ms2_threshold))

            if tmp_df.shape[0] == 1:

                obs_score_df = obs_score_df.append(tmp_df)
                obs_fp_lst.append(mz)
                fp_obs_lst.append(1)
            elif tmp_df.shape[0] > 1:
                tmp_df = tmp_df.sort_values(by='i', ascending=False)
                obs_score_df = obs_score_df.append(tmp_df.head(1))
                obs_fp_lst.append(mz)
                fp_obs_lst.append(1)
            else:
                tmp_df = pd.DataFrame(data={'mz': [0.0], 'i': [0.0]})
                obs_score_df = obs_score_df.append(tmp_df)
                missed_fp_lst.append(mz)
                fp_obs_lst.append(0)

        obs_lst = obs_score_df['mz'].tolist()

        # fingerprint_score = 100 * (1 - spatial.distance.cosine(np.array(obs_lst), np.array(fingerprint_lst)))
        print(fp_lib_lst)
        print(fp_obs_lst)
        fp_sim_score = 100 * (1 - spatial.distance.cosine(np.array(obs_lst), np.array(fingerprint_lst)))
        fp_sim_score = round(fp_sim_score, 1)
        print('fp_sim_score', fp_sim_score)
        fingerprint_score = 100 * (1 - spatial.distance.cosine(np.array(fp_obs_lst), np.array(fp_lib_lst)))
        fingerprint_score = round(fingerprint_score, 1)
        print('fingerprint_score', fingerprint_score)

        fp_info_dct = {'fingerprint_score': fingerprint_score, 'obs_score_df': obs_score_df,
                       'obs_mz': obs_fp_lst, 'missed_mz': missed_fp_lst}

        return fp_info_dct

    def get_specific_peaks(self, mz_lib, ms2_df, ms2_precision=50e-6, ms2_threshold=10,
                           ms2_hginfo_threshold=0.02, vendor='waters'):

        _target_frag_df = pd.DataFrame()
        _target_nl_df = pd.DataFrame()
        _other_frag_df = pd.DataFrame()
        _other_nl_df = pd.DataFrame()

        ms2_max_i = ms2_df['i'].max()
        ms2_hginfo_abs_i = ms2_max_i * ms2_hginfo_threshold
        ms2_threshold = max(ms2_threshold, ms2_hginfo_abs_i)

        for _i, _frag_se in self.target_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']
            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0
                _delta = _frag_mz * ms2_precision + seg_shift
                _frag_mz_low = _frag_mz - _delta
                _frag_mz_high = _frag_mz + _delta
            else:
                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _target_frag_df = _target_frag_df.append(_frag_df.head(1))

        for _i, _frag_se in self.other_frag_df.iterrows():

            _frag_mz = _frag_se['EXACTMASS']
            _frag_class = _frag_se['CLASS']
            _frag_label = _frag_se['LABEL']

            if vendor == 'thermo':
                if _frag_mz < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = _frag_mz * ms2_precision + seg_shift
                _frag_mz_low = _frag_mz - _delta
                _frag_mz_high = _frag_mz + _delta

            else:
                _frag_mz_low = _frag_mz - _frag_mz * ms2_precision
                _frag_mz_high = _frag_mz + _frag_mz * ms2_precision
            _frag_mz_query_code = '%f <= mz <= %f and i > %f' % (_frag_mz_low, _frag_mz_high, ms2_threshold)

            _frag_df = ms2_df.query(_frag_mz_query_code)

            if _frag_df.shape[0] > 0:
                _frag_df = _frag_df.sort_values(by='i', ascending=False)
                _frag_df.loc[:, 'CLASS'] = _frag_class
                _frag_df.loc[:, 'LABEL'] = _frag_label
                _frag_df.loc[:, _frag_label] = 100 * _frag_df['i'] / ms2_max_i
                _other_frag_df = _other_frag_df.append(_frag_df.head(1))

        for _i, _nl_se in self.target_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']
            _nl_label = _nl_se['LABEL']

            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = (mz_lib - _nl_mz) * ms2_precision * ms2_precision + seg_shift
                _nl_mz_low = mz_lib - _nl_mz - _delta
                _nl_mz_high = mz_lib - _nl_mz + _delta

            else:
                _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
                _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision

            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _nl_df.loc[:, 'LABEL'] = _nl_label
                _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
                _target_nl_df = _target_nl_df.append(_nl_df.head(1))

        for _i, _nl_se in self.other_nl_df.iterrows():

            _nl_mz = _nl_se['EXACTMASS']
            _nl_class = _nl_se['CLASS']
            _nl_label = _nl_se['LABEL']

            if vendor == 'thermo':
                if mz_lib < 450.0:
                    seg_shift = 0.0
                elif 450.0 <= mz_lib < 600.0:
                    seg_shift = 0.1
                elif 600.0 <= mz_lib < 750.0:
                    seg_shift = 0.2
                elif 750.0 <= mz_lib < 900.0:
                    seg_shift = 0.3
                elif 900.0 <= mz_lib < 1200.0:
                    seg_shift = 0.4
                elif 1200.0 <= mz_lib:
                    seg_shift = 0.5
                else:
                    seg_shift = 0.0

                _delta = (mz_lib - _nl_mz) * ms2_precision + seg_shift
                _nl_mz_low = mz_lib - _nl_mz - _delta
                _nl_mz_high = mz_lib - _nl_mz + _delta

            else:
                _nl_mz_low = mz_lib - _nl_mz - _nl_mz * ms2_precision
                _nl_mz_high = mz_lib - _nl_mz + _nl_mz * ms2_precision

            _nl_mz_query_code = '%f <= mz <= %f and i > %f' % (_nl_mz_low, _nl_mz_high, ms2_threshold)

            _nl_df = ms2_df.query(_nl_mz_query_code)

            if _nl_df.shape[0] > 0:
                _nl_df = _nl_df.sort_values(by='i', ascending=False)
                _nl_df.loc[:, 'CLASS'] = _nl_class
                _nl_df.loc[:, 'LABEL'] = _nl_label
                _nl_df.loc[:, _nl_label] = 100 * _nl_df['i'] / ms2_max_i
                _other_nl_df = _other_nl_df.append(_nl_df.head(1))

        specific_ion_dct = {}
        if _target_frag_df.shape[0] > 0:
            specific_ion_dct['TARGET_FRAG'] = _target_frag_df
        if _target_nl_df.shape[0] > 0:
            specific_ion_dct['TARGET_NL'] = _target_nl_df
        if _other_frag_df.shape[0] > 0:
            specific_ion_dct['OTHER_FRAG'] = _other_frag_df
        if _other_nl_df.shape[0] > 0:
            specific_ion_dct['OTHER_NL'] = _other_nl_df

        return specific_ion_dct

    @staticmethod
    def get_snr_score(ident_info_dct, specific_ion_dct, msp_obs_df, fp_obs_df, amplify_factor=3.5767):

        signal_df = pd.DataFrame()
        tmp_s_df = ident_info_dct['MATCHED_FA_INFO']
        tmp_s_df = pd.DataFrame(tmp_s_df, columns=['mz', 'i'])
        signal_df = signal_df.append(tmp_s_df)
        try:
            tmp_s_df = ident_info_dct['MATCHED_LYSO_INFO']
            tmp_s_df = pd.DataFrame(tmp_s_df, columns=['mz', 'i'])
            signal_df = signal_df.append(tmp_s_df)
        except KeyError:
            pass
        try:
            tmp_s_df = specific_ion_dct['TARGET_FRAG']
            tmp_s_df = pd.DataFrame(tmp_s_df, columns=['mz', 'i'])
            signal_df = signal_df.append(tmp_s_df)
        except KeyError:
            pass
        try:
            tmp_s_df = specific_ion_dct['TARGET_NL']
            tmp_s_df = pd.DataFrame(tmp_s_df, columns=['mz', 'i'])
            signal_df = signal_df.append(tmp_s_df)
        except KeyError:
            pass
        try:
            tmp_s_df = pd.DataFrame(msp_obs_df, columns=['mz', 'i'])
            signal_df = signal_df.append(tmp_s_df)
        except KeyError:
            pass
        # remove fingerprint for signal consideration
        # try:
        #     tmp_s_df = pd.DataFrame(fp_obs_df, columns=['mz', 'i'])
        #     signal_df = signal_df.append(tmp_s_df)
        # except KeyError:
        #     pass

        noise_df = pd.DataFrame()
        try:
            tmp_n_df = ident_info_dct['FA_INFO']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass
        try:
            tmp_n_df = ident_info_dct['LYSO_INFO']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass
        try:
            tmp_n_df = ident_info_dct['LYSO_W_INFO']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass
        try:
            tmp_n_df = ident_info_dct['OTHER_SIGNALS_INFO']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass
        try:
            tmp_n_df = specific_ion_dct['OTHER_FRAG']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass
        try:
            tmp_n_df = specific_ion_dct['OTHER_NL']
            tmp_n_df = pd.DataFrame(tmp_n_df, columns=['mz', 'i'])
            noise_df = noise_df.append(tmp_n_df)
        except KeyError:
            pass

        signal_df = signal_df.drop_duplicates()
        signal_sum_i = sum(signal_df['i'].tolist())

        print('signal_df')
        print(signal_df)

        noise_min_i = 0

        if noise_df.shape[0] > 0:

            noise_df = noise_df.drop_duplicates()
            # remove mz if they are identified
            noise_df = noise_df[~noise_df.isin(signal_df.to_dict('l'))]
            if noise_df.shape[0] > 0:
                noise_df = noise_df.loc[:, ['mz', 'i']]
                noise_df = noise_df.dropna(how='any')
                if noise_df.shape[0] > 0:
                    noise_sum_i = sum(noise_df['i'].tolist())
                    noise_min_i = noise_df['i'].min()
                    if noise_sum_i == 0:
                        noise_sum_i = 1
                else:
                    noise_sum_i = 1
            else:
                noise_sum_i = 1
            print('noise_df')
            print(noise_df)
        else:
            noise_sum_i = 1

        # use the SNR equation SNR = 20 * log10(signal/noise)
        # snr_score = 20 * math.log10((signal_sum_i / noise_sum_i))
        # set s/n == 25 --> SNR_SCORE = 100
        # default 3.5767 = 100 / (20 * math.log10(25)) --> 3.5767

        sn_ratio = signal_sum_i / noise_sum_i

        snr_score = amplify_factor * 20 * math.log10(sn_ratio)

        if snr_score < 0:
            snr_score = 0.0
        elif snr_score > 100.0:
            snr_score = 100.0
        else:
            snr_score = round(snr_score, 1)

        print('sn_ratio ==>', sn_ratio)
        print('SNR SCORE ==>', snr_score)

        snr_i_info = {'min_signal': signal_df['i'].min(), 'min_noise': noise_min_i}

        return snr_score, sn_ratio, noise_df, snr_i_info


if __name__ == '__main__':
    pass
