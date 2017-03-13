# -*- coding: utf-8 -*-
#
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPtiger.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#


from __future__ import division

import json
import math
import os
import time
from multiprocessing import Pool

import pandas as pd

from LibHunter.SpectraExtractor import extract_mzml
from LibHunter.SpectraExtractor import get_spectra
from LibHunter.SpectraExtractor import get_xic_all
from LibHunter.ScoreGenerator import ScoreGenerator
from LibHunter.PanelPlotter import plot_spectra
from LibHunter.IsotopeHunter import IsotopeHunter
from LibHunter.LogPageCreator import LogPageCreator
from LibHunter.PrecursorHunter import PrecursorHunter
from LibHunter.ScoreHunter import get_lpp_info


def huntlipids(param_dct):
    """

    :param param_dct:
    :return:
    """

    start_time = time.clock()

    usr_core_num = 4

    usr_lipid_type = param_dct['lipid_type']
    charge_mode = param_dct['charge_mode']
    usr_vendor = param_dct['vendor']

    usr_xlsx = param_dct['lpp_sum_info_path_str']
    # usr_sdf = param_dct['sdf_path_str']
    # usr_msp = param_dct['msp_path_str']
    usr_mzml = param_dct['mzml_path_str']
    output_folder = param_dct['img_output_folder_str']
    output_sum_xlsx = param_dct['xlsx_output_path_str']

    fa_list_cfg = param_dct['fa_sum_path_str']
    key_frag_cfg = param_dct['lipid_specific_cfg']
    score_cfg = param_dct['score_cfg']

    usr_rt_range = [param_dct['rt_start'], param_dct['rt_end']]
    mz_start = param_dct['mz_start']
    mz_end = param_dct['mz_end']
    usr_pr_mz_range = [mz_start, mz_end]
    usr_dda_top = param_dct['dda_top']
    usr_ms1_threshold = param_dct['ms_th']
    usr_ms2_threshold = param_dct['ms2_th']
    usr_ms2_hg_threshold = param_dct['hg_th']
    usr_ms1_precision = param_dct['ms_ppm'] * 1e-6
    usr_ms2_precision = param_dct['ms2_ppm'] * 1e-6
    usr_ms2_hg_precision = param_dct['hg_ppm'] * 1e-6
    usr_pr_window = param_dct['pr_window']

    usr_rankscore_filter = param_dct['score_filter']
    usr_isotope_score_filter = param_dct['isotope_score_filter']
    usr_ms2_info_th = param_dct['ms2_infopeak_threshold']
    usr_ms2_hginfo_th = param_dct['ms2_hginfopeak_threshold']
    usr_rank_mode = param_dct['rank_score']
    usr_fast_isotope = param_dct['fast_isotope']

    # use the SNR equation SNR = 20 * log10(signal/noise)
    # snr_score = 20 * math.log10((signal_sum_i / noise_sum_i))
    # set s/n == 20 --> SNR_SCORE = 100
    # default 6.5051 = 100 / (20 * math.log10(20)) --> 6.5051
    usr_max_sn_ratio = param_dct['sn_ratio']
    if usr_max_sn_ratio == 20 or usr_max_sn_ratio == 0:
        usr_amp_factor = 6.5051
    else:
        usr_amp_factor = 100 / (20 * math.log10(usr_max_sn_ratio))

    if usr_fast_isotope is True:
        isotope_score_mode = '(Fast mode)'
    else:
        isotope_score_mode = ''

    usr_ms1_ppm = int(param_dct['ms_ppm'])

    hunter_start_time_str = param_dct['hunter_start_time']
    isotope_hunter = IsotopeHunter()

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
    # usr_fa_def_df.loc[:, 'C'] = usr_fa_def_df['C'].astype(int)
    # usr_fa_def_df.loc[:, 'DB'] = usr_fa_def_df['DB'].astype(int)

    usr_weight_df = pd.read_excel(score_cfg, index_col='Type')
    usr_weight_df = usr_weight_df.loc[:, 'Weight']

    usr_key_frag_df = pd.read_excel(key_frag_cfg)
    usr_key_frag_df = usr_key_frag_df.query('EXACTMASS > 0')

    # get the information from the following columns and leave the rewark back
    usr_key_frag_df = usr_key_frag_df[['CLASS', 'TYPE', 'EXACTMASS', 'PR_CHARGE', 'LABEL', 'CHARGE_MODE']]

    print('=== ==> --> Start to parse mzML')
    # extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                    ms1_threshold=usr_ms1_threshold, ms2_threshold=usr_ms2_threshold,
                                                    ms1_precision=usr_ms1_precision, ms2_precision=usr_ms2_precision,
                                                    vendor=usr_vendor
                                                    )

    ms1_obs_pr_df = pr_hunter.get_matched_pr(usr_scan_info_df, usr_spectra_pl, core_num=usr_core_num)

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

    ms1_xic_mz_lst = ms1_obs_pr_df['MS1_XIC_mz'].tolist()
    ms1_xic_mz_lst = set(ms1_xic_mz_lst)

    print('=== ==> --> Start to extract XIC')
    sub_len = int(math.ceil(len(ms1_xic_mz_lst) / usr_core_num))
    core_key_list = map(None, *(iter(ms1_xic_mz_lst),) * sub_len)
    print(core_key_list)
    # Start multiprocessing
    print('!!!!!! Start multiprocessing ==> ==> ==> Number of Cores: %i' % usr_core_num)
    xic_dct = {}
    parallel_pool = Pool()
    xic_results_lst = []
    for core_list in core_key_list:
        core_list = filter(lambda x: x is not None, core_list)
        xic_result = parallel_pool.apply_async(get_xic_all, args=(core_list, usr_mzml, usr_rt_range,
                                                                  usr_ms1_precision, usr_ms2_precision,
                                                                  usr_vendor,))
        xic_results_lst.append(xic_result)

    parallel_pool.close()
    parallel_pool.join()

    for xic_result in xic_results_lst:
        sub_xic_dct = xic_result.get()
        if len(sub_xic_dct.keys()) > 0:
            if len(xic_dct.keys()) > 0:
                xic_dct = dict(xic_dct, **sub_xic_dct)
            else:
                xic_dct = sub_xic_dct.copy()

    print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

    target_ident_lst = []
    ident_page_idx = 1
    checked_info_df.sort_values(by=['Lib_mz', 'scan_time', 'MS2_PR_mz'],
                                ascending=[True, True, True], inplace=True)

    print('=== ==> --> Start to Hunt for LPPs !!')
    checked_info_groups = checked_info_df.groupby(['MS2_PR_mz', 'Lib_mz', 'Formula', 'scan_time'])
    all_group_key_lst = checked_info_groups.groups.keys()
    sub_len = int(math.ceil(len(all_group_key_lst) / usr_core_num))
    core_key_list = map(None, *(iter(all_group_key_lst),) * sub_len)
    print(core_key_list)
    # Start multiprocessing
    print('!!!!!! Start multiprocessing ==> ==> ==> Number of Cores: %i' % usr_core_num)
    parallel_pool = Pool()
    lpp_info_results_lst = []
    for core_list in core_key_list:
        core_list = filter(lambda x: x is not None, core_list)
        lpp_info_result = parallel_pool.apply_async(get_lpp_info, args=(param_dct, checked_info_df,
                                                                        checked_info_groups, core_list,
                                                                        usr_fa_def_df, usr_weight_df, usr_key_frag_df,
                                                                        usr_scan_info_df, ms1_xic_mz_lst,
                                                                        usr_spectra_pl, xic_dct, target_ident_lst))
        lpp_info_results_lst.append(lpp_info_result)

    parallel_pool.close()
    parallel_pool.join()

    for lpp_info_result in lpp_info_results_lst:
        tmp_lpp_info_df = lpp_info_result.get()

        if tmp_lpp_info_df.shape[0] > 0:

            output_df = output_df.append(tmp_lpp_info_df)

    ident_page_idx += 1

    # log_pager.add_info(img_name_core, ident_page_idx, _tmp_output_df)

    print('=== ==> --> Generate the output table')
    if output_df.shape[0] > 0:
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

        output_header_lst = ['Bulk_identification', 'Proposed_structures', 'Formula_neutral', 'Formula_ion',
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
        output_df.to_excel(output_sum_xlsx, index=False)
        print(output_sum_xlsx)
        print('=== ==> --> saved >>> >>> >>>')

    log_pager.close_page()

    tot_run_time = time.clock() - start_time

    tot_run_time = '%.2f Sec\n' % tot_run_time

    print('>>> >>> >>> FINISHED in %s sec <<< <<< <<<' % tot_run_time)

    return tot_run_time


if __name__ == '__main__':
    pass
