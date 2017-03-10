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

import pandas as pd

from LibHunter.SpectraExtractor import extract_mzml
from LibHunter.SpectraExtractor import get_spectra
from LibHunter.SpectraExtractor import get_xic_all
from LibHunter.ScoreGenerator import ScoreGenerator
from LibHunter.PanelPlotter import plot_spectra
from LibHunter.ScoreFilter import check_peaks
from LibHunter.IsotopeHunter import IsotopeHunter
# from LibHunter.AbbrElemCalc import BulkAbbrFormula
from LibHunter.LogPageCreator import LogPageCreator
from LibHunter.PrecursorHunter import PrecursorHunter


def huntlipids(param_dct):
    """

    :param param_dct:
    :return:
    """

    start_time = time.clock()

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

    # hunter_folder = param_dct['hunter_folder']

    if usr_rank_mode is True:
        score_mode = 'Rank mode'
    else:
        score_mode = 'Relative intensity mode'

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

    pr_hunter = PrecursorHunter(lpp_info_df, usr_mzml, param_dct)
    # abbr2formula = BulkAbbrFormula()

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

    score_calc = ScoreGenerator(usr_fa_def_df, usr_weight_df, usr_key_frag_df, usr_lipid_type, ion_charge=charge_mode)

    print('=== ==> --> Start to parse mzML')
    # extract all spectra from mzML to pandas DataFrame
    usr_scan_info_df, usr_spectra_pl = extract_mzml(usr_mzml, usr_rt_range, dda_top=usr_dda_top,
                                                    ms1_threshold=usr_ms1_threshold, ms2_threshold=usr_ms2_threshold,
                                                    ms1_precision=usr_ms1_precision, ms2_precision=usr_ms2_precision,
                                                    vendor=usr_vendor
                                                    )

    ms1_obs_pr_df = pr_hunter.get_matched_pr(usr_scan_info_df, usr_spectra_pl)

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

    ms1_obs_mz_lst = ms1_obs_pr_df['MS1_obs_mz'].tolist()
    ms1_obs_mz_lst = set(ms1_obs_mz_lst)

    ms1_xic_mz_lst = ms1_obs_pr_df['MS1_XIC_mz'].tolist()
    ms1_xic_mz_lst = set(ms1_xic_mz_lst)

    print('=== ==> --> Start to extract XIC')
    try:
        xic_dct = get_xic_all(ms1_obs_pr_df, usr_mzml, usr_rt_range, ms1_precision=usr_ms1_precision,
                              msn_precision=usr_ms2_precision, vendor=usr_vendor)

    except KeyError:
        return u'Nothing found! Check mzML vendor settings!'

    print('=== ==> --> Number of XIC extracted: %i' % len(xic_dct.keys()))

    # plot_info_dct = {}
    # ms1_pr_mz_lst = []
    target_ident_lst = []
    ident_page_idx = 1
    print('checked_info_df')
    print(checked_info_df.shape)

    # get spectra of one ABBR and plot
    for _n, _subgroup_df in checked_info_df.groupby(['MS2_PR_mz', 'Lib_mz', 'Formula', 'scan_time', 'Abbreviation']):
        _row_se = _subgroup_df.iloc[0, :].squeeze()
        _usr_ms2_pr_mz = _row_se['MS2_PR_mz']
        # _usr_ms1_obs_mz = _row_se['MS1_obs_mz']
        _usr_ms2_rt = _row_se['scan_time']
        # _usr_formula_charged = _row_se['Formula']
        _usr_charge = _row_se['Ion']
        _usr_ms2_dda_rank = _row_se['DDA_rank']
        _usr_ms2_scan_id = _row_se['scan_number']
        _usr_mz_lib = _row_se['Lib_mz']
        _usr_abbr_bulk = _row_se['Abbreviation']
        _tmp_chk_df = usr_scan_info_df.query('MS2_PR_mz == %.6f and DDA_rank == %i and scan_number == %i'
                                             % (_usr_ms2_pr_mz, _usr_ms2_dda_rank, _usr_ms2_scan_id))

        _usr_formula = _row_se['FORMULA_NEUTRAL']
        _usr_formula_charged = _row_se['Formula']

        if _tmp_chk_df.shape[0] == 1:
            print('>>> >>> >>> Processing:', _tmp_chk_df.head())
            print('>>> >>> >>> >>> MS2 PR m/z %f' % _usr_ms2_pr_mz)
            usr_spec_info_dct = get_spectra(_usr_ms2_pr_mz, _usr_mz_lib, _usr_ms2_dda_rank, _usr_ms2_scan_id,
                                            ms1_xic_mz_lst, usr_scan_info_df, usr_spectra_pl,
                                            dda_top=usr_dda_top, ms1_precision=usr_ms1_precision, vendor=usr_vendor
                                            )
            _ms1_pr_i = usr_spec_info_dct['ms1_i']
            _ms1_pr_mz = usr_spec_info_dct['ms1_mz']
            _ms1_df = usr_spec_info_dct['ms1_df']
            _ms2_df = usr_spec_info_dct['ms2_df']

            if _ms1_pr_mz > 0.0 and _ms1_df.shape[0] > 0 and _ms2_df.shape[0] > 0 and _ms1_pr_i > 0.0:

                print('>>> >>> >>> >>> Best PR on MS1: %f' % _ms1_pr_mz)

                isotope_score_info_dct = isotope_hunter.get_isotope_score(_ms1_pr_mz, _ms1_pr_i,
                                                                          _usr_formula_charged, _ms1_df,
                                                                          isotope_number=2,
                                                                          only_c=usr_fast_isotope,
                                                                          score_filter=usr_isotope_score_filter)

                isotope_score = isotope_score_info_dct['isotope_score']

                print('isotope_score: %f' % isotope_score)
                if isotope_score >= usr_isotope_score_filter:
                    print('>>> isotope_check PASSED! >>> >>> >>>')
                    print('>>> >>> >>> >>> Entry Info >>> >>> >>> >>> ')
                    _row_se.set_value('MS1_obs_mz', _ms1_pr_mz)
                    _exact_ppm = 1e6 * (_ms1_pr_mz - _usr_mz_lib) / _usr_mz_lib
                    _row_se.set_value('ppm', _exact_ppm)
                    _row_se.set_value('abs_ppm', abs(_exact_ppm))
                    _msp_df = pd.read_json(_row_se['MSP_JSON'], orient='index')
                    print(_row_se)
                    match_info_dct, matched_checker = score_calc.get_rankscore(_usr_abbr_bulk, charge_mode, _ms1_pr_mz,
                                                                               _ms2_df,
                                                                               ms2_precision=usr_ms2_precision,
                                                                               ms2_threshold=usr_ms2_threshold,
                                                                               ms2_infopeak_threshold=usr_ms2_info_th)
                    rank_score = match_info_dct['Rank_score']
                    matched_fa_df = match_info_dct['MATCHED_FA_INFO']

                    if matched_checker > 0 and rank_score > usr_rankscore_filter and matched_fa_df.shape[0] > 0:
                        print('Rank_score: %.f --> passed' % rank_score)

                        fingerprint_lst = json.loads(_row_se['FINGERPRINT'])
                        fp_info_dct = score_calc.get_fingerprint_score(fingerprint_lst, _ms2_df,
                                                                       ms2_precision=usr_ms2_precision,
                                                                       ms2_threshold=usr_ms2_threshold,
                                                                       ms2_infopeak_threshold=usr_ms2_info_th)
                        _fp_score = fp_info_dct['fingerprint_score']
                        _obs_fp_df = fp_info_dct['obs_score_df']
                        obs_fp_lst = fp_info_dct['obs_mz']
                        missed_fp_lst = fp_info_dct['missed_mz']

                        if _fp_score > 0:

                            _cosine_score, _msp_df, _obs_msp_df = score_calc.get_cosine_score(_msp_df, _ms2_df,
                                                                                              ms2_precision=
                                                                                              usr_ms2_precision,
                                                                                              ms2_threshold=
                                                                                              usr_ms2_threshold,
                                                                                              ms2_infopeak_threshold=
                                                                                              usr_ms2_info_th
                                                                                              )

                            print('==> --> Cosine similarity score: %f' % _cosine_score)

                            if _cosine_score > 0:

                                specific_check_dct = score_calc.get_specific_peaks(_usr_mz_lib, _ms2_df,
                                                                                   ms2_precision=
                                                                                   usr_ms2_hg_precision,
                                                                                   ms2_threshold=
                                                                                   usr_ms2_hg_threshold,
                                                                                   ms2_hginfo_threshold=
                                                                                   usr_ms2_hginfo_th,
                                                                                   vendor=usr_vendor
                                                                                   )

                                # format abbr. for file names
                                _save_abbr_bulk = _usr_abbr_bulk
                                _save_abbr_bulk = _save_abbr_bulk.replace(r'(', r'[')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r')', r']')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r'<', r'[')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r'>', r']')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r':', r'-')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r'@', r'-')
                                _save_abbr_bulk = _save_abbr_bulk.replace('\\', r'_')
                                _save_abbr_bulk = _save_abbr_bulk.replace(r'/', r'_')

                                img_name_core = ('\%.4f_rt%.3f_DDAtop%.0f_scan%.0f_%s.png'
                                                 % (_usr_ms2_pr_mz, _usr_ms2_rt, _usr_ms2_dda_rank,
                                                    _usr_ms2_scan_id, _save_abbr_bulk)
                                                 )

                                img_name = (output_folder + r'\LPPtiger_Results_Figures_%s' % hunter_start_time_str
                                            + img_name_core)

                                snr_score, sn_ratio, noise_df, snr_i_info = score_calc.get_snr_score(match_info_dct,
                                                                                                     specific_check_dct,
                                                                                                     _obs_msp_df,
                                                                                                     _obs_fp_df,
                                                                                                     amplify_factor=
                                                                                                     usr_amp_factor)

                                if sn_ratio >= 1.0:

                                    overall_score = sum([rank_score, _cosine_score, _fp_score,
                                                         snr_score, isotope_score]) / 5
                                    overall_score = round(overall_score, 1)

                                    if overall_score >= 40:
                                        match_info_dct['Cosine_score'] = _cosine_score
                                        match_info_dct['Fingerprint'] = _fp_score
                                        match_info_dct['Isotope_score'] = isotope_score
                                        match_info_dct['SNR_score'] = snr_score
                                        match_info_dct['Overall_score'] = overall_score

                                        isotope_checker, isotope_score = plot_spectra(_row_se, xic_dct,
                                                                                      match_info_dct,
                                                                                      usr_spec_info_dct,
                                                                                      specific_check_dct,
                                                                                      isotope_score_info_dct,
                                                                                      _usr_formula_charged, _usr_charge,
                                                                                      save_img_as=img_name,
                                                                                      ms1_precision=usr_ms1_precision,
                                                                                      msp_info=_msp_df,
                                                                                      obs_fp=obs_fp_lst,
                                                                                      missed_fp=missed_fp_lst,
                                                                                      snr_i_info=snr_i_info
                                                                                      )

                                        print('==> check for output -->')

                                        if _ms1_pr_i > 0 and isotope_checker == 0 \
                                                and isotope_score > usr_isotope_score_filter:

                                            if 'OTHER_FRAG' in specific_check_dct.keys():
                                                other_frag_df = specific_check_dct['OTHER_FRAG']
                                                other_frag_count = other_frag_df.shape[0]
                                            else:
                                                other_frag_count = 0
                                            if 'OTHER_NL' in specific_check_dct.keys():
                                                other_nl_df = specific_check_dct['OTHER_NL']
                                                other_nl_count = other_nl_df.shape[0]
                                            else:
                                                other_nl_count = 0
                                            if 'TARGET_FRAG' in specific_check_dct.keys():
                                                target_frag_df = specific_check_dct['TARGET_FRAG']
                                                target_frag_count = target_frag_df.shape[0]
                                                target_frag_col_lst = target_frag_df.columns.tolist()
                                                for _frag_abbr in target_frag_col_lst:
                                                    if _frag_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                                                        for _i, _f_se in target_frag_df.iterrows():
                                                            if _f_se['LABEL'] == _frag_abbr:
                                                                match_info_dct[_frag_abbr] = _f_se[_frag_abbr]
                                                                if _frag_abbr not in target_ident_lst:
                                                                    target_ident_lst.append(_frag_abbr)
                                            else:
                                                target_frag_count = 0
                                            if 'TARGET_NL' in specific_check_dct.keys():
                                                target_nl_df = specific_check_dct['TARGET_NL']
                                                target_nl_count = target_nl_df.shape[0]
                                                target_nl_col_lst = target_nl_df.columns.tolist()
                                                for _nl_abbr in target_nl_col_lst:
                                                    if _nl_abbr not in ['mz', 'i', 'LABEL', 'CLASS']:
                                                        for _i, _n_se in target_nl_df.iterrows():
                                                            if _n_se['LABEL'] == _nl_abbr:
                                                                match_info_dct[_nl_abbr] = _n_se[_nl_abbr]
                                                                if _nl_abbr not in target_ident_lst:
                                                                    target_ident_lst.append(_nl_abbr)
                                            else:
                                                target_nl_count = 0

                                            match_info_dct['Proposed_structures'] = _usr_abbr_bulk
                                            match_info_dct['Bulk_identification'] = _usr_abbr_bulk
                                            match_info_dct['Formula_neutral'] = _usr_formula
                                            match_info_dct['Formula_ion'] = _usr_formula_charged
                                            match_info_dct['Charge'] = _usr_charge
                                            match_info_dct['MS1_obs_mz'] = _ms1_pr_mz
                                            match_info_dct['MS1_obs_i'] = '%.2e' % float(_ms1_pr_i)
                                            match_info_dct['Lib_mz'] = _usr_mz_lib
                                            match_info_dct['MS2_scan_time'] = _usr_ms2_rt
                                            match_info_dct['DDA#'] = _usr_ms2_dda_rank
                                            match_info_dct['MS2_PR_mz'] = _usr_ms2_pr_mz
                                            match_info_dct['Scan#'] = _usr_ms2_scan_id
                                            match_info_dct['#Specific_peaks'] = target_frag_count + target_nl_count
                                            match_info_dct['#Contaminated_peaks'] = other_frag_count + other_nl_count
                                            match_info_dct['ppm'] = _exact_ppm
                                            match_info_dct['SN_ratio'] = '%.1f' % sn_ratio

                                            try:
                                                del match_info_dct['MATCH_INFO']
                                                del match_info_dct['OTHER_SIGNALS_INFO']
                                                del match_info_dct['MATCHED_FA_INFO']
                                                del match_info_dct['MATCHED_LYSO_INFO']
                                            except KeyError:
                                                pass

                                            print('match_info_dct')
                                            print(match_info_dct)

                                            _tmp_output_df = pd.DataFrame(data=match_info_dct,
                                                                          index=[ident_page_idx])

                                            output_df = output_df.append(_tmp_output_df)

                                            log_pager.add_info(img_name_core, ident_page_idx, _tmp_output_df)
                                            ident_page_idx += 1

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
