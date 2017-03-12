# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPsmi.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from __future__ import division
from __future__ import print_function

import gc

import pandas as pd
import numpy as np

from ParallelFunc import ppm_calc_para, ppm_window_para, pr_window_calc_para


class PrecursorHunter(object):
    def __init__(self, lpp_info_df, param_dct):
        self.lpp_info_df = lpp_info_df
        self.param_dct = param_dct
        gc.disable()

    def get_matched_pr(self, scan_info_df, spectra_pl):

        print('Start match!!!!')

        pr_window = self.param_dct['pr_window']
        ms1_th = self.param_dct['ms_th']
        ms1_ppm = self.param_dct['ms_ppm']

        pl_class = self.param_dct['lipid_type']
        charge_mode = self.param_dct['charge_mode']

        ms1_obs_pr_df = pd.DataFrame()

        if pl_class == 'PC' and charge_mode == '[M+HCOO]-':

            pc_fa_df = self.lpp_info_df[self.lpp_info_df['[M+HCOO]-_MZ'] > 0]
            pc_h_df = self.lpp_info_df[self.lpp_info_df['[M-H]-_MZ'] > 0]

            print('PC [M+HCOO]- LPP: ', pc_fa_df.shape[0])
            print('PC [M-H]- LPP: ', pc_h_df.shape[0])
            pc_fa_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(pc_fa_df['[M+HCOO]-_MZ'].values, -1 * pr_window)
            pc_fa_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(pc_fa_df['[M+HCOO]-_MZ'].values, pr_window)
            pc_fa_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(pc_fa_df['[M+HCOO]-_MZ'].values, -1 * ms1_ppm)
            pc_fa_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(pc_fa_df['[M+HCOO]-_MZ'].values, ms1_ppm)
            pc_fa_df.loc[:, 'Formula'] = pc_fa_df['[M+HCOO]-_FORMULA'].str.strip('-')
            pc_fa_df.loc[:, 'Ion'] = '[M+HCOO]-'
            pc_fa_df.loc[:, 'Lib_mz'] = pc_fa_df.loc[:, '[M+HCOO]-_MZ']
            if pc_h_df.shape[0] > 0:
                pc_h_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(pc_h_df['[M-H]-_MZ'].values, -1 * pr_window)
                pc_h_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(pc_h_df['[M-H]-_MZ'].values, pr_window)
                pc_h_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(pc_h_df['[M-H]-_MZ'].values, -1 * ms1_ppm)
                pc_h_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(pc_h_df['[M-H]-_MZ'].values, ms1_ppm)
                pc_h_df.loc[:, 'Formula'] = pc_h_df['[M-H]-_FORMULA'].str.strip('-')
                pc_h_df.loc[:, 'Ion'] = '[M-H]-'
                pc_h_df.loc[:, 'Lib_mz'] = pc_h_df.loc[:, '[M-H]-_MZ']
                self.lpp_info_df = pc_fa_df.copy()
                self.lpp_info_df = self.lpp_info_df.append(pc_h_df)
            else:
                self.lpp_info_df = pc_fa_df.copy()

            self.lpp_info_df = self.lpp_info_df.sort_values(by='PR_MZ_LOW')

        else:
            self.lpp_info_df.loc[:, 'PR_MZ_LOW'] = pr_window_calc_para(self.lpp_info_df['[M-H]-_MZ'].values,
                                                                       -1 * pr_window)
            self.lpp_info_df.loc[:, 'PR_MZ_HIGH'] = pr_window_calc_para(self.lpp_info_df['[M-H]-_MZ'].values,
                                                                        pr_window)
            self.lpp_info_df.loc[:, 'MS1_MZ_LOW'] = ppm_window_para(self.lpp_info_df['[M-H]-_MZ'].values, -1 * ms1_ppm)
            self.lpp_info_df.loc[:, 'MS1_MZ_HIGH'] = ppm_window_para(self.lpp_info_df['[M-H]-_MZ'].values, ms1_ppm)
            self.lpp_info_df.loc[:, 'Formula'] = self.lpp_info_df.loc[:, '[M-H]-_FORMULA'].str.strip('-')
            self.lpp_info_df.loc[:, 'Ion'] = '[M-H]-'
            self.lpp_info_df.loc[:, 'Lib_mz'] = self.lpp_info_df.loc[:, '[M-H]-_MZ']

        for _n, _subgroup_df in self.lpp_info_df.groupby(['Lib_mz', 'Formula']):
            _samemz_se = _subgroup_df.iloc[0, :].squeeze()

            _pr_code = '%f<= MS2_PR_mz <= %f' % (_samemz_se['PR_MZ_LOW'], _samemz_se['PR_MZ_HIGH'])

            _tmp_scan_info_df = scan_info_df.query(_pr_code)

            if _tmp_scan_info_df.shape[0] > 0:

                for idx, row in _tmp_scan_info_df.iterrows():

                    ms2_pr = row['MS2_PR_mz']
                    ms2_dda_idx = row['dda_event_idx']
                    ms2_dda_rank = row['DDA_rank']
                    ms2_spec_idx = row['spec_index']
                    ms2_scan_time = row['scan_time']
                    ms2_scan_number = row['scan_number']

                    _ms1_query_code = 'dda_event_idx == %i and DDA_rank == 0' % ms2_dda_idx
                    tmp_ms1_info_df = scan_info_df.query(_ms1_query_code)

                    if tmp_ms1_info_df.shape[0] > 0:

                        ms1_spec_idx = tmp_ms1_info_df['spec_index'].tolist()[0]

                        if ms1_spec_idx in spectra_pl.items:
                            ms1_df = spectra_pl[ms1_spec_idx]
                            pr_ms1_df = ms1_df.query('i > %f and %f <= mz <= %f' % (ms1_th, _samemz_se['MS1_MZ_LOW'],
                                                                                    _samemz_se['MS1_MZ_HIGH']))
                            if pr_ms1_df.shape[0] > 0:

                                # find the best ms1 pr with min abs_ppm
                                # pr_ms1_df.loc[:, 'ppm'] = (1e6 * (pr_ms1_df.loc[:, 'mz'] - _samemz_se['Lib_mz']) /
                                #                            _samemz_se['Lib_mz'])
                                pr_ms1_df.loc[:, 'ppm'] = ppm_calc_para(pr_ms1_df['mz'].values, _samemz_se['Lib_mz'])
                                pr_ms1_df.loc[:, 'abs_ppm'] = abs(pr_ms1_df.loc[:, 'ppm'])
                                pr_info_df = pr_ms1_df.query('abs_ppm <= %i' % ms1_ppm)
                                if pr_info_df.shape[0] > 0:
                                    pr_info_df = pr_info_df.sort_values(by='i', ascending=False).reset_index(drop=True)
                                    _ms1_mz = pr_info_df['mz'].tolist()[0]
                                    _ppm = pr_info_df['ppm'].tolist()[0]
                                    _subgroup_df.loc[:, 'MS1_obs_mz'] = _ms1_mz
                                    _subgroup_df.loc[:, 'dda_event_idx'] = ms2_dda_idx
                                    _subgroup_df.loc[:, 'spec_index'] = ms2_spec_idx
                                    _subgroup_df.loc[:, 'scan_time'] = ms2_scan_time
                                    _subgroup_df.loc[:, 'DDA_rank'] = ms2_dda_rank
                                    _subgroup_df.loc[:, 'scan_number'] = ms2_scan_number
                                    _subgroup_df.loc[:, 'MS2_PR_mz'] = ms2_pr
                                    _subgroup_df.loc[:, 'MS1_XIC_mz'] = round(_ms1_mz, 4)
                                    _subgroup_df.loc[:, 'ppm'] = _ppm
                                    _subgroup_df.loc[:, 'abs_ppm'] = abs(_ppm)
                                    # print(_subgroup_df[['Abbreviation', 'spec_index', 'MS2_PR_mz', 'MS1_obs_mz']])
                                    ms1_obs_pr_df = ms1_obs_pr_df.append(_subgroup_df)

        # print('ms1_obs_pr_df.shape', ms1_obs_pr_df.shape)
        if ms1_obs_pr_df.shape[0] > 0:
            ms1_obs_pr_df = ms1_obs_pr_df.sort_values(by=['Lib_mz', 'abs_ppm'], ascending=[True, True])
            ms1_obs_pr_df = ms1_obs_pr_df.reset_index(drop=True)
            return ms1_obs_pr_df

        else:
            return '!! NO suitable precursor --> Check settings!!'
