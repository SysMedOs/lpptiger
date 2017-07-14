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

from __future__ import division

import json
import math
import os
import time
import sys
from multiprocessing import Pool

import pandas as pd

from LibHunter.SpectraExtractor import extract_mzml
from LibHunter.SpectraExtractor import get_spectra
from LibHunter.SpectraExtractor import get_xic_from_pl
from LibHunter.SpectraExtractor import get_spec_info
from LibHunter.ScoreGenerator import ScoreGenerator
from LibHunter.PanelPlotter import plot_spectra
from LibHunter.IsotopeHunter import IsotopeHunter
from LibHunter.LogPageCreator import LogPageCreator
from LibHunter.PrecursorHunter import PrecursorHunter
from LibHunter.ScoreHunter import get_lpp_info


def huntlipids(param_dct):
    """
    hunter_param_dct = {'hunter_folder': self.theolpp_cwd, 'hunter_start_time': start_time_str,
                            'vendor': usr_vendor, 'Experiment_mode': usr_exp_mode,
                            'lipid_type': _pl_class, 'charge_mode': _pl_charge,
                            'lpp_sum_info_path_str': lpp_sum_info_path_str, 'fa_sum_path_str': fa_sum_path_str,
                            'mzml_path_str': mzml_path_str, 'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str,
                            'rt_start': rt_start, 'rt_end': rt_end, 'mz_start': mz_start, 'mz_end': mz_end,
                            'ms_th': ms_th, 'ms_ppm': ms_ppm, 'pr_window': pr_window, 'dda_top': dda_top,
                            'ms2_th': ms2_th, 'ms2_ppm': ms2_ppm, 'ms2_infopeak_threshold': ms2_info_threshold,
                            'hg_th': hg_th, 'hg_ppm': hg_ppm, 'ms2_hginfopeak_threshold': hgms2_info_threshold,
                            'score_filter': overall_score_filter,
                            'lipid_specific_cfg': lipid_specific_cfg, 'score_cfg': score_cfg, 'sn_ratio': sn_ratio,
                            'isotope_score_filter': isotope_score_filter, 'rank_score_filter': rank_score_filter,
                            'msp_score_filter': msp_score_filter, 'fp_score_filter': fp_score_filter,
                            'snr_score_filter': snr_score_filter,
                            'parallization_mode': parallization_mode, 'core_number': core_num, 'max_ram': max_ram,
                            'img_type': img_typ, 'img_dpi': img_dpi, 'fast_isotope': fast_isotope,
                            }
    :param param_dct:
    :return:
    """

    start_time = time.clock()

    hunter_start_time_str = param_dct['hunter_start_time']
    usr_vendor = param_dct['vendor']

    usr_lipid_type = param_dct['lipid_type']
    usr_xlsx = param_dct['lpp_sum_info_path_str']
    usr_mzml = param_dct['mzml_path_str']
    output_folder = param_dct['img_output_folder_str']
    output_sum_xlsx = param_dct['xlsx_output_path_str']
    fa_list_cfg = param_dct['fa_sum_path_str']
    key_frag_cfg = param_dct['lipid_specific_cfg']
    score_cfg = param_dct['score_cfg']
    usr_rt_range = [param_dct['rt_start'], param_dct['rt_end']]
    mz_start = param_dct['mz_start']
    mz_end = param_dct['mz_end']
    usr_dda_top = param_dct['dda_top']
    usr_ms1_threshold = param_dct['ms_th']
    usr_ms1_max = param_dct['ms_max']
    usr_ms2_threshold = param_dct['ms2_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_precision = param_dct['ms2_ppm'] * 1e-6

    # parameters from settings tab
    usr_core_num = param_dct['core_number']
    usr_max_ram = param_dct['max_ram']

    lpp_info_df = pd.read_excel(usr_xlsx)

    # cut lib info to the user defined m/z range
    if usr_lipid_type == 'PC':
        tmp_lpp_info_df = lpp_info_df[(mz_start <= lpp_info_df['[M-H]-_MZ']) & (lpp_info_df['[M-H]-_MZ'] <= mz_end)]
        tmp_fa_lpp_info_df = lpp_info_df[(mz_start <= lpp_info_df['[M+HCOO]-_MZ']) &
                                         (lpp_info_df['[M+HCOO]-_MZ'] <= mz_end)]
        lpp_info_df = tmp_lpp_info_df.append(tmp_fa_lpp_info_df).copy()
    else:
        lpp_info_df = lpp_info_df[(mz_start <= lpp_info_df['[M-H]-_MZ']) & (lpp_info_df['[M-H]-_MZ'] <= mz_end)]

    pr_hunter = PrecursorHunter(lpp_info_df, param_dct)

    # keep stay in current working directory
    current_path = os.getcwd()
    if os.path.isdir(output_folder):
        os.chdir(output_folder)
        if os.path.isdir('LPPtiger_Results_Figures_%s' % hunter_start_time_str):
            pass
        else:
            os.mkdir('LPPtiger_Results_Figures_%s' % hunter_start_time_str)
    os.chdir(current_path)

    log_pager = LogPageCreator(output_folder, hunter_start_time_str, param_dct)

    output_df = pd.DataFrame()

    print('=== ==> --> Start to process')
    print('=== ==> --> Phospholipid class: %s' % usr_lipid_type)

    # generate the indicator table

    usr_fa_def_df = pd.read_excel(fa_list_cfg)

    usr_weight_df = pd.read_excel(score_cfg, index_col='Type')
    usr_weight_df = usr_weight_df.loc[:, 'Weight']

    usr_key_frag_df = pd.read_excel(key_frag_cfg)
    usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')

    # get the information from the following columns and leave the rewark back
    usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]

    print('=== ==> --> Start to parse mzML')
    # extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl, ms1_xic_df = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                                ms1_threshold=usr_ms1_threshold,
                                                                ms2_threshold=usr_ms2_threshold,
                                                                ms1_precision=usr_ms1_precision,
                                                                ms2_precision=usr_ms2_precision,
                                                                vendor=usr_vendor, ms1_max=usr_ms1_max
                                                                )

    print('MS1_XIC_df.shape', ms1_xic_df.shape)

    ms1_obs_pr_df, sub_pl_group_lst = pr_hunter.get_matched_pr(usr_scan_info_df, usr_spectra_pl, ms1_max=usr_ms1_max,
                                                               core_num=usr_core_num, max_ram=usr_max_ram)

    if isinstance(ms1_obs_pr_df, str):
        return '!! NO suitable precursor --> Check settings!!\n'

    print('=== ==> --> ms1 precursor matched')

    # remove bad precursors
    checked_info_df = pd.DataFrame()
    for _idx, _check_scan_se in usr_scan_info_df.iterrows():
        _dda_rank = _check_scan_se['DDA_rank']
        _scan_id = _check_scan_se['scan_number']
        _tmp_usr_df = ms1_obs_pr_df.query('DDA_rank == %f and scan_number == %f' % (_dda_rank, _scan_id))
        checked_info_df = checked_info_df.append(_tmp_usr_df)

    checked_info_df.sort_values(by='MS2_PR_mz')

    ms1_xic_mz_lst = ms1_obs_pr_df['MS1_XIC_mz'].tolist()
    ms1_xic_mz_lst = sorted(set(ms1_xic_mz_lst))
    print('ms1_xic_mz_lst', len(ms1_xic_mz_lst))
    print(ms1_xic_mz_lst)

    print('=== ==> --> Start to extract XIC')
    if len(ms1_xic_mz_lst) >= 20 * usr_core_num:
        sub_len = int(math.ceil(len(ms1_xic_mz_lst) / usr_core_num))
        core_key_list = map(None, *(iter(ms1_xic_mz_lst),) * sub_len)
    else:
        core_key_list = [ms1_xic_mz_lst]
    # print(core_key_list)
    # Start multiprocessing
    print('!!!!!! Start multiprocessing ==> ==> ==> Number of Cores: %i' % usr_core_num)
    xic_dct = {}

    # if usr_core_num >= 8:
    #     xic_core_num = 8
    # else:
    #     xic_core_num = usr_core_num
    # parallel_pool = Pool(xic_core_num)
    if 1 < usr_core_num < len(core_key_list):
        parallel_pool = Pool(usr_core_num)
        xic_results_lst = []
        core_worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = filter(lambda x: x is not None, core_list)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                print(core_list)
                # xic_result = parallel_pool.apply_async(get_xic_all, args=(core_list, usr_mzml, usr_rt_range,
                #                                                           usr_ms1_precision, usr_ms2_precision,
                #                                                           usr_vendor,))
                xic_result = parallel_pool.apply_async(get_xic_from_pl, args=(core_list, ms1_xic_df, 500))
                core_worker_count += 1
                xic_results_lst.append(xic_result)

        parallel_pool.close()
        parallel_pool.join()

        for xic_result in xic_results_lst:
            try:
                sub_xic_dct = xic_result.get()
                if len(sub_xic_dct.keys()) > 0:
                    xic_dct = dict(xic_dct, **sub_xic_dct)
            except (KeyError, SystemError, ValueError):
                pass
    else:
        print('Using single core mode...')
        core_worker_count = 1
        for core_list in core_key_list:
            if isinstance(core_list, tuple) or isinstance(core_list, list):
                if None in core_list:
                    core_list = filter(lambda x: x is not None, core_list)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                print(core_list)
                sub_xic_dct = get_xic_from_pl(core_list, ms1_xic_df, 500)
                core_worker_count += 1
                if len(sub_xic_dct.keys()) > 0:
                    xic_dct = dict(xic_dct, **sub_xic_dct)

    print('xic_dct', len(xic_dct.keys()))
    print(xic_dct.keys())

    if len(xic_dct.keys()) == 0:
        print('No precursor for XIC found')
        return '!! NO suitable precursor --> Check settings!!\n'
    else:
        print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

    target_ident_lst = []
    # ident_page_idx = 1
    checked_info_df.sort_values(by=['Lib_mz', 'scan_time', 'MS2_PR_mz'],
                                ascending=[True, True, True], inplace=True)

    print('=== ==> --> Start to Hunt for LPPs !!')
    checked_info_groups = checked_info_df.groupby(['Lib_mz', 'MS2_PR_mz', 'Formula', 'scan_time', 'Ion'])
    lpp_all_group_key_lst = checked_info_groups.groups.keys()
    # lpp_all_group_key_lst = sorted(lpp_all_group_key_lst, key=lambda x: x[0])

    spec_sub_len = int(math.ceil(len(lpp_all_group_key_lst) / usr_core_num))
    spec_sub_key_lst = map(None, *(iter(lpp_all_group_key_lst),) * spec_sub_len)

    lpp_spec_info_dct = {}

    if usr_core_num > 1:
        parallel_pool = Pool(usr_core_num)
        spec_results_lst = []
        core_worker_count = 1
        for _sub_lst in spec_sub_key_lst:
            if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                if None in _sub_lst:
                    _sub_lst = filter(lambda x: x is not None, _sub_lst)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                spec_result = parallel_pool.apply_async(get_spec_info, args=(_sub_lst, checked_info_groups,
                                                                             usr_scan_info_df))
                core_worker_count += 1
                spec_results_lst.append(spec_result)

        parallel_pool.close()
        parallel_pool.join()

        for spec_result in spec_results_lst:
            try:
                sub_spec_dct = spec_result.get()
                if len(sub_spec_dct.keys()) > 0:
                    lpp_spec_info_dct = dict(lpp_spec_info_dct, **sub_spec_dct)
            except (KeyError, SystemError, ValueError):
                print('ValueError: must supply a tuple to get_group with multiple grouping keys')
    else:
        print('Using single core mode...')
        core_worker_count = 1
        for _sub_lst in spec_sub_key_lst:
            if isinstance(_sub_lst, tuple) or isinstance(_sub_lst, list):
                if None in _sub_lst:
                    _sub_lst = filter(lambda x: x is not None, _sub_lst)
                else:
                    pass
                print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                sub_spec_dct = get_spec_info(_sub_lst, checked_info_groups, usr_scan_info_df)
                core_worker_count += 1
                if len(sub_spec_dct.keys()) > 0:
                    lpp_spec_info_dct = dict(lpp_spec_info_dct, **sub_spec_dct)

    print('lpp_spec_info_dct', len(lpp_spec_info_dct.keys()))

    # Single process ONLY. usr_spectra_pl is too big in RAM --> RAM leaking during copy
    lpp_spec_dct = {}
    spec_info_key_lst = lpp_spec_info_dct.keys()
    for _spec_group_key in spec_info_key_lst:
        _spec_info_dct = lpp_spec_info_dct[_spec_group_key]
        _usr_ms2_pr_mz = _spec_info_dct['MS2_PR_mz']
        _usr_ms2_dda_rank = _spec_info_dct['DDA_rank']
        _usr_ms2_scan_id = _spec_info_dct['scan_number']
        _usr_mz_lib = _spec_info_dct['Lib_mz']
        usr_spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_dda_rank, _usr_ms2_scan_id,
                                        ms1_xic_mz_lst, usr_scan_info_df, usr_spectra_pl,
                                        dda_top=usr_dda_top, ms1_precision=usr_ms1_precision, vendor=usr_vendor
                                        )
        lpp_spec_dct[_spec_group_key] = usr_spec_info_dct

    found_spec_key_lst = lpp_spec_dct.keys()
    found_spec_key_lst = sorted(found_spec_key_lst, key=lambda x: x[0])
    spec_key_num = len(found_spec_key_lst)
    lpp_part_key_lst = []
    if spec_key_num > (usr_core_num * 40):
        lpp_part_len = int(math.ceil(spec_key_num / 8))
        lpp_part_lst = map(None, *(iter(found_spec_key_lst),) * lpp_part_len)
        for part_lst in lpp_part_lst:
            if None in part_lst:
                part_lst = filter(lambda x: x is not None, part_lst)
            lpp_sub_len = int(math.ceil(len(part_lst) / usr_core_num))
            lpp_sub_key_lst = map(None, *(iter(part_lst),) * lpp_sub_len)
            lpp_part_key_lst.append(lpp_sub_key_lst)

    else:
        lpp_sub_len = int(math.ceil(spec_key_num / usr_core_num))
        lpp_sub_key_lst = map(None, *(iter(found_spec_key_lst),) * lpp_sub_len)
        lpp_part_key_lst.append(lpp_sub_key_lst)

    part_tot = len(lpp_part_key_lst)
    part_counter = 1

    for lpp_sub_key_lst in lpp_part_key_lst:

        if part_tot == 1:
            print('>>> Start multiprocessing ==> Max Number of Cores: %i' % usr_core_num)
        else:
            print('>>> Start multiprocessing ==> Part %i / %i --> Max Number of Cores: %i' %
                  (part_counter, part_tot, usr_core_num))
        part_counter += 1
        # Start multiprocessing
        if usr_core_num > 1:
            parallel_pool = Pool(usr_core_num)
            lpp_info_results_lst = []
            core_worker_count = 1
            for lpp_sub_lst in lpp_sub_key_lst:
                if isinstance(lpp_sub_lst, tuple) or isinstance(lpp_sub_lst, list):
                    if None in lpp_sub_lst:
                        lpp_sub_lst = filter(lambda x: x is not None, lpp_sub_lst)
                    else:
                        pass
                    lpp_sub_dct = {k: lpp_spec_dct[k] for k in lpp_sub_lst}
                    print('>>> >>> Core #%i ==> ...... processing ......' % core_worker_count)
                    lpp_info_result = parallel_pool.apply_async(get_lpp_info, args=(param_dct, checked_info_df,
                                                                                    checked_info_groups, lpp_sub_lst,
                                                                                    usr_fa_def_df, usr_weight_df,
                                                                                    usr_key_frag_df,
                                                                                    usr_scan_info_df, ms1_xic_mz_lst,
                                                                                    lpp_sub_dct, xic_dct, target_ident_lst))
                    # ('>>> >>> Get lpp_info_result of this worker-->', <class 'multiprocessing.pool.ApplyResult'>)
                    lpp_info_results_lst.append(lpp_info_result)
                    core_worker_count += 1

            parallel_pool.close()
            parallel_pool.join()

            for lpp_info_result in lpp_info_results_lst:
                try:
                    tmp_lpp_info_df = lpp_info_result.get()
                except (KeyError, SystemError, ValueError):
                    tmp_lpp_info_df = 'error'
                    print('!!error!!--> This segment receive no LPP identified.')
                if isinstance(tmp_lpp_info_df, str):
                    pass
                else:
                    if isinstance(tmp_lpp_info_df, pd.DataFrame):
                        if tmp_lpp_info_df.shape[0] > 0:
                            output_df = output_df.append(tmp_lpp_info_df)
        else:
            print('Using single core mode...')
            core_worker_count = 1
            for lpp_sub_lst in lpp_sub_key_lst:

                if isinstance(lpp_sub_lst, tuple) or isinstance(lpp_sub_lst, list):
                    if None in lpp_sub_lst:
                        lpp_sub_lst = filter(lambda x: x is not None, lpp_sub_lst)
                    else:
                        pass
                    lpp_sub_dct = {k: lpp_spec_dct[k] for k in lpp_sub_lst}
                    print('>>> >>> Part %i Subset #%i ==> ...... processing ......' % (part_counter, core_worker_count))
                    tmp_lpp_info_df = get_lpp_info(param_dct, checked_info_df, checked_info_groups, lpp_sub_lst,
                                                   usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_scan_info_df,
                                                   ms1_xic_mz_lst, lpp_sub_dct, xic_dct, target_ident_lst)
                    core_worker_count += 1
                    if isinstance(tmp_lpp_info_df, str):
                        pass
                    else:
                        if tmp_lpp_info_df.shape[0] > 0:
                            output_df = output_df.append(tmp_lpp_info_df)

    print('=== ==> --> Generate the output table')
    if output_df.shape[0] > 0:
        try:
            output_df = output_df.sort_values(by=['Lib_mz', 'Proposed_structures', 'MS2_scan_time', 'Overall_score'])
        except KeyError:
            pass
        output_df.reset_index(drop=True, inplace=True)
        output_df.index += 1
        # print('output_df')
        # print(output_df.head(5))
        # print(output_df.columns.tolist())
        output_df.drop_duplicates(keep='first', inplace=True)
        log_pager.add_all_info(output_df)
        output_header_lst = output_df.columns.tolist()
        for _i_check in ['i_sn1', 'i_sn2', 'i_[M-H]-sn1', 'i_[M-H]-sn2', 'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O']:
            if _i_check not in output_header_lst:
                output_df[_i_check] = 0.0

        output_round_dct = {r'MS1_obs_mz': 4, r'Lib_mz': 4, 'ppm': 2, 'MS2_scan_time': 3,
                            'i_sn1': 2, 'i_sn2': 2, 'i_[M-H]-sn1': 2, 'i_[M-H]-sn2': 2,
                            'i_[M-H]-sn1-H2O': 2, 'i_[M-H]-sn2-H2O': 2
                            }
        # add intensities of target peaks to round list
        if len(target_ident_lst) > 0:
            for _t in target_ident_lst:
                output_round_dct[_t] = 2
        output_df = output_df.round(output_round_dct)

        output_df.rename(columns={'#Contaminated_peaks': '#Unspecific_peaks'}, inplace=True)

        output_header_lst = ['Proposed_structures', 'Formula_neutral', 'Formula_ion',
                             'Charge', 'Lib_mz', 'ppm', 'SN_ratio', 'Overall_score',
                             'Rank_score', 'Cosine_score', 'Fingerprint', 'SNR_score', 'Isotope_score',
                             'MS1_obs_mz', 'MS1_obs_i', r'MS2_PR_mz', 'MS2_scan_time',
                             'DDA#', 'Scan#', 'i_sn1', 'i_sn2',
                             'i_[M-H]-sn1', 'i_[M-H]-sn2', 'i_[M-H]-sn1-H2O', 'i_[M-H]-sn2-H2O', '#Specific_peaks']
        output_header_lst += target_ident_lst
        output_header_lst += ['#Unspecific_peaks']

        output_df = output_df[output_header_lst]
        output_df = output_df.sort_values(by=['MS1_obs_mz', 'MS2_scan_time', 'Rank_score'],
                                          ascending=[True, True, False])
        output_df = output_df.reset_index(drop=True)
        output_df.index += 1
        try:
            output_df.to_excel(output_sum_xlsx, index=False)
            print(output_sum_xlsx)
        except IOError:
            output_df.to_excel('%s-%i%s' % (output_sum_xlsx[:-5], int(time.time()), '.xlsx'), index=False)
            print(output_sum_xlsx)
        print('=== ==> --> saved >>> >>> >>>')

    log_pager.close_page()

    tot_run_time = time.clock() - start_time

    tot_run_time = '%.2f Sec\n' % tot_run_time

    print('>>> >>> >>> FINISHED in %s sec <<< <<< <<<' % tot_run_time)

    return tot_run_time


if __name__ == '__main__':
    pass
