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
import os
import json

from PySide import QtCore, QtGui
import matplotlib

# matplotlib.use('Qt4Agg')
matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import pandas as pd


def plot_spectra(abbr, mz_se, xic_dct, ident_info_dct, spec_info_dct, specific_check_dct, isotope_score_info_dct,
                 formula_charged, charge, save_img_as=None, img_type='png', dpi=300,
                 ms1_precision=50e-6, msp_info=pd.DataFrame(), obs_fp=[], missed_fp=[], snr_i_info={}):
    if msp_info.shape[0] > 0:
        plot_msp = 1
    else:
        plot_msp = 0

    if plot_msp == 1:
        msp_abs_min = abs(max(msp_info['rev_abs_i'].tolist()))
    else:
        msp_abs_min = 0

    ms2_pr_mz = mz_se['MS2_PR_mz']
    ms1_obs = mz_se['MS1_obs_mz']
    ms1_xic_mz = mz_se['MS1_XIC_mz']
    lib_mz = mz_se['Lib_mz']
    func_id = mz_se['DDA_rank']
    ms1_pr_ppm = mz_se['ppm']

    isotope_score = isotope_score_info_dct['isotope_score']
    isotope_checker_dct = isotope_score_info_dct['isotope_checker_dct']
    m2_score = isotope_score_info_dct['m2_score']
    m2_checker_dct = isotope_score_info_dct['m2_checker_dct']
    deconv_lst = isotope_score_info_dct['deconv_lst']
    if len(deconv_lst) == 4:
        pass
    else:
        deconv_lst = [0, 0, 0, 0]

    ms1_delta = lib_mz * ms1_precision

    ms1_pr_i = spec_info_dct['ms1_i']
    ms1_pr_mz = spec_info_dct['ms1_mz']
    ms1_rt = spec_info_dct['ms1_rt']
    ms2_rt = spec_info_dct['ms2_rt']
    ms1_df = spec_info_dct['ms1_df']
    ms2_df = spec_info_dct['ms2_df']

    ms1_df = ms1_df.sort_values(by='mz', ascending='True')
    dash_i = [ms1_df['i'].max()]

    ms_zoom_query_str = ' %.2f < mz < %.2f' % (ms1_obs - 1.5, ms1_obs + 3.55)
    ms_zoom_df = ms1_df.query(ms_zoom_query_str)
    try:
        ms_zoom_bp_i = max(ms_zoom_df['i'].tolist())
    except ValueError:
        ms_zoom_bp_i = 0

    xic_df = xic_dct[ms1_xic_mz]

    xic_rt_lst = xic_df['rt'].tolist()
    xic_i_lst = xic_df['i'].tolist()

    if ms_zoom_bp_i > 0 and len(xic_rt_lst) > 0 and len(xic_i_lst) > 0:
        # cut lower peaks to accelerate plotting time
        m1_dct = isotope_checker_dct[1]
        m1_theo_mz = m1_dct['theo_mz']
        m1_theo_i = m1_dct['theo_i']
        m1_obs_mz = m1_dct['obs_mz']
        m1_obs_i = m1_dct['obs_i']
        print('spectra shape', ms1_df.shape, ms2_df.shape)
        if ms1_df['i'].max() >= 10000 and ms1_df.shape[0] >= 500:
            ms1_min = ms1_df['i'].min()
            ms1_max = ms1_df['i'].max()
            ms1_top1000_i = sorted(ms1_df['i'].tolist(), reverse=True)[499]
            ms1_plot_th = min(m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
            ms1_plot_th = max(ms1_plot_th, ms1_top1000_i)
            print(m1_obs_i, 3 * ms1_min, ms1_max * 0.01, 1000, ms1_top1000_i)
            ms1_df = ms1_df.query('i >= %f' % ms1_plot_th)
            print('Plot full MS1 with abs intensity filter > %f' % ms1_plot_th)
        if ms2_df['i'].max() >= 1000 and ms2_df.shape[0] >= 500:
            ms2_min = ms2_df['i'].min()
            ms2_max = ms2_df['i'].max()

            ms2_top1000_i = sorted(ms2_df['i'].tolist(), reverse=True)[499]
            ms2_min_lst = [3 * ms2_min, ms2_max * 0.01, 10, msp_abs_min, ms2_top1000_i]
            ms2_plot_th = max(min(ms2_min_lst), ms2_top1000_i)
            try:
                ms2_sig_i_lst = []
                for _sn in snr_i_info.keys():
                    ms2_sig_i_lst.append(snr_i_info[_sn])

                ms2_sig_th = min(ms2_sig_i_lst)
                if ms2_sig_th < ms2_plot_th:
                    ms2_plot_th = ms2_sig_th
            except KeyError:
                pass
            print(ms2_min_lst)
            ms2_plot_th -= 1
            if ms2_plot_th > 0:
                ms2_df = ms2_df.query('i >= %f' % ms2_plot_th)
                print('Plot full MS/MS with abs intensity filter > %f' % ms2_plot_th)

        _msms_low_df = ms2_df.query('mz <= 400')
        _msms_high_df = ms2_df.query('mz > 400')
        _msms_high_df = _msms_high_df.query('mz < %.4f' % (ms2_pr_mz + 1))

        print ('Start looking for MS2 PR m/z %f @ MS1 best PR m/z %f with lib m/z %f'
               % (ms2_pr_mz, ms1_obs, lib_mz))

        # Generate A4 image in landscape
        fig, pic_array = plt.subplots(nrows=3, ncols=2, figsize=(11.692, 8.267), sharex=False,
                                      sharey=False)
        # Make better spacing between subplots
        plt.tight_layout()
        xic_pic = pic_array[0, 0]
        msms_pic = pic_array[0, 1]
        ms_pic = pic_array[1, 0]
        msms_low_pic = pic_array[1, 1]
        ms_zoom_pic = pic_array[2, 0]
        msms_high_pic = pic_array[2, 1]

        xic_pic.tick_params(axis='both', which='major', labelsize=10)
        msms_pic.tick_params(axis='both', which='major', labelsize=10)
        ms_pic.tick_params(axis='both', which='major', labelsize=10)
        msms_low_pic.tick_params(axis='both', which='major', labelsize=10)
        ms_zoom_pic.tick_params(axis='both', which='major', labelsize=10)
        msms_high_pic.tick_params(axis='both', which='major', labelsize=10)

        # ms spectrum start
        # ms_pic.stem(ms1_df['mz'].tolist(), ms1_df['i'].tolist(), 'grey', markerfmt=' ')
        ms_pic.plot(ms1_df['mz'].tolist(), ms1_df['i'].tolist(), 'grey')

        markerline, stemlines, baseline = ms_pic.stem([ms1_pr_mz], dash_i, markerfmt=' ')
        plt.setp(stemlines, color='#00ccff', linewidth=5, alpha=0.3)
        markerline, stemlines, baseline = ms_pic.stem([ms1_pr_mz], [ms1_pr_i], markerfmt='D')
        plt.setp(stemlines, color='magenta')
        plt.setp(markerline, markerfacecolor='magenta', markersize=4, markeredgewidth=0)

        ms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        ms_pic.set_ylabel("Intensity", fontsize=10)
        ms_pic.set_ylim([0, max(ms1_df['i'].tolist()) * 1.3])

        # add annotation
        _ms_pkl_top_df = ms1_df.sort_values(by='i', ascending=False).head(10)
        _ms_pkl_top_peak_list = zip(_ms_pkl_top_df['mz'].tolist(), _ms_pkl_top_df['i'].tolist())
        for _ms_pkl_top_peak in _ms_pkl_top_peak_list:
            _ms_pkl_top_peak_str = '%.4f' % _ms_pkl_top_peak[0]
            _ms_pkl_top_peak_y = _ms_pkl_top_peak[1]
            ms_pic.text(_ms_pkl_top_peak[0], _ms_pkl_top_peak_y, _ms_pkl_top_peak_str, fontsize=6)

        m0_theo_base_box = patches.Rectangle((lib_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[0],
                                             facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m0_theo_base_box)
        m0_theo_box = patches.Rectangle((lib_mz - ms1_delta, deconv_lst[0]), 2 * ms1_delta, ms1_pr_i - deconv_lst[0],
                                        facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m0_theo_box)

        # isotope region | if any peak in M-1.0034

        m_pre_theo_box = patches.Rectangle((lib_mz - 1.0034 - ms1_delta, 0), 2 * ms1_delta, ms1_pr_i,
                                           facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none')
        ms_zoom_pic.add_patch(m_pre_theo_box)

        ms_zoom_offset_i = ms_zoom_bp_i * 0.1

        ms_zoom_pic.set_xlim([ms1_pr_mz - 1.5, ms1_pr_mz + 3.55])
        ms_zoom_pic.set_ylim([0, ms_zoom_bp_i * 1.45])
        ms_zoom_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), fontsize=10)
        ms_zoom_pic.ticklabel_format(axis='x', useOffset=False, fontsize=10)
        ms_zoom_pic.set_xlabel('m/z', fontsize=10, labelpad=-1)
        ms_zoom_pic.set_ylabel('Intensity', fontsize=10)

        # ms_zoom_pic.stem(ms_zoom_df['mz'].tolist(), ms_zoom_df['i'].tolist(), 'grey', markerfmt=' ', zorder=1)
        ms_zoom_df = ms_zoom_df.sort_values(by='mz', ascending='True')
        ms_zoom_pic.plot(ms_zoom_df['mz'].tolist(), ms_zoom_df['i'].tolist(), 'grey', zorder=1)
        markerline, stemlines, baseline = ms_zoom_pic.stem([ms1_pr_mz], [ms1_pr_i],
                                                           'magenta', markerfmt='D', zorder=20)
        plt.setp(markerline, markerfacecolor='magenta', markeredgecolor='none', markeredgewidth=0,
                 markersize=6, alpha=0.8)
        ms_zoom_pic.text(ms1_pr_mz + 0.06, ms1_pr_i, '%.4f' % float(ms1_pr_mz),
                         color='magenta', fontsize=6
                         )
        markerline, stemlines, baseline = ms_zoom_pic.stem([lib_mz], [ms1_pr_i], '--', markerfmt='o', zorder=21)
        plt.setp(markerline, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0)
        plt.setp(stemlines, color='#ff6600')
        ms_zoom_pic.text(lib_mz - 0.15, ms1_pr_i + ms_zoom_offset_i, '[M+0]', color='#ff6600', fontsize=6)
        ms_zoom_pic.text(lib_mz - 0.71, ms1_pr_i, 'Calc: %.4f' % lib_mz, color='#ff6600', fontsize=6)

        # isotope region | highlight the 1st isotope

        # theo range box
        m1_theo_base_box = patches.Rectangle((m1_theo_mz - ms1_delta, 0),
                                             2 * ms1_delta, deconv_lst[1],
                                             facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m1_theo_base_box)
        # m1_theo_box = patches.Rectangle((m1_theo_mz - ms1_delta, deconv_lst[1]), 2 * ms1_delta, m1_theo_i - deconv_lst[1],
        #                                 facecolor=(0.1, 1.0, 1.0, 0.3), edgecolor='none', zorder=7)
        m1_theo_box = patches.Rectangle((m1_theo_mz - ms1_delta, deconv_lst[1]), 2 * ms1_delta, m1_theo_i - deconv_lst[1],
                                        facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
        ms_zoom_pic.add_patch(m1_theo_box)

        markerline, stemlines, baseline = ms_zoom_pic.stem([m1_theo_mz], [m1_theo_i], '--', markerfmt='o', zorder=22)
        plt.setp(stemlines, color='#ff6600')
        plt.setp(markerline, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0)
        ms_zoom_pic.text(m1_theo_mz - 0.15, m1_theo_i + ms_zoom_offset_i, '[M+1]', color='#ff6600', fontsize=6)
        ms_zoom_pic.text(m1_theo_mz - 0.71, m1_theo_i, 'Calc: %.4f' % m1_theo_mz,
                         color='#ff6600', fontsize=6)
        ms_zoom_pic.text(m1_obs_mz + 0.04, m1_obs_i, '%.4f' % m1_obs_mz, color='magenta', fontsize=6)

        opt_box_lst = []

        # isotope region | highlight the 2nd isotope
        if 2 in isotope_checker_dct.keys():
            m2_dct = isotope_checker_dct[2]
            m2_theo_mz = m2_dct['theo_mz']
            m2_theo_i = m2_dct['theo_i']
            m2_obs_mz = m2_dct['obs_mz']
            m2_obs_i = m2_dct['obs_i']
            # m2_theo_r = m2_dct['theo_ratio']
            # # m2_obs_r = m2_dct['obs_ratio']
            # m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_i,
            #                                 facecolor=(0.2, 1.0, 1.0, 0.3), edgecolor='none', zorder=9)
            m2_theo_box = patches.Rectangle((m2_theo_mz - ms1_delta, 0), 2 * ms1_delta, m2_theo_i,
                                            facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
            ms_zoom_pic.add_patch(m2_theo_box)
            opt_box_lst.append(ms_zoom_pic)
            markerline, stemlines, baseline = ms_zoom_pic.stem([m2_theo_mz], [m2_theo_i], '--', markerfmt='o', zorder=23)
            plt.setp(stemlines, color='#ff6600', alpha=0.8)
            plt.setp(markerline, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0, alpha=0.9)
            ms_zoom_pic.text(m2_theo_mz - 0.15, m2_theo_i + ms_zoom_offset_i, '[M+2]', color='#ff6600', fontsize=6)
            plt.setp(markerline, markerfacecolor='#ff6600', markersize=6, markeredgewidth=0, alpha=0.9)
            ms_zoom_pic.text(m2_theo_mz - 0.71, m2_theo_i, 'Calc: %.4f' % m2_theo_mz,
                             color='#ff6600', fontsize=6)
            ms_zoom_pic.text(m2_obs_mz + 0.04, m2_obs_i, '%.4f' % m2_obs_mz, color='magenta', fontsize=6)

        if len(m2_checker_dct.keys()) > 0:
            for _mh2 in m2_checker_dct.keys():
                mh2_dct = m2_checker_dct[_mh2]
                mh2_theo_mz = mh2_dct['theo_mz']
                mh2_theo_i = mh2_dct['theo_i']
                mh2_obs_mz = mh2_dct['obs_mz']
                mh2_obs_i = mh2_dct['obs_i']
                decon_idx = _mh2 + 2
                # mh2_theo_r = mh2_dct['theo_ratio']
                # mh2_obs_r = mh2_dct['obs_ratio']
                mh2_theo_base_box = patches.Rectangle((mh2_theo_mz - ms1_delta, 0), 2 * ms1_delta, deconv_lst[decon_idx],
                                                      facecolor=(0, 0.8, 1.0, 0.6), edgecolor='none', zorder=1)
                ms_zoom_pic.add_patch(mh2_theo_base_box)
                # opt_box_lst.append(mh2_theo_base_box)
                # mh2_theo_box = patches.Rectangle((mh2_theo_mz - ms1_delta, deconv_lst[decon_idx]),
                #                                  2 * ms1_delta, mh2_theo_i - deconv_lst[decon_idx],
                #                                  facecolor=(1.0, 0.0, 0.0, 0.4), edgecolor='none', zorder=12)
                mh2_theo_box = patches.Rectangle((mh2_theo_mz - ms1_delta, deconv_lst[decon_idx]),
                                                 2 * ms1_delta, mh2_theo_i - deconv_lst[decon_idx],
                                                 facecolor=(1.0, 0.0, 0.0, 0.6), edgecolor='none', zorder=1)
                ms_zoom_pic.add_patch(mh2_theo_box)
                markerline, stemlines, baseline = ms_zoom_pic.stem([mh2_theo_mz], [mh2_theo_i], '--',
                                                                   markerfmt='o', zorder=24)
                plt.setp(stemlines, color='red', alpha=0.8)
                plt.setp(markerline, markerfacecolor='red', markersize=6, markeredgewidth=0, alpha=0.9)
                if _mh2 == 0:
                    _mh2_name = ''
                else:
                    _mh2_name = '+%i' % _mh2
                ms_zoom_pic.text(mh2_theo_mz - 0.2 - 0.05 * _mh2, mh2_theo_i + ms_zoom_offset_i,
                                 '[M+H2%s]' % _mh2_name, color='red', fontsize=6)
                ms_zoom_pic.text(mh2_theo_mz - 0.71, mh2_theo_i, 'Calc: %.4f' % mh2_theo_mz,
                                 color='red', fontsize=6)
                ms_zoom_pic.text(mh2_obs_mz + 0.04, mh2_obs_i, '%.4f' % mh2_obs_mz, color='red', fontsize=6)

            # plot the M+H2 isotope score

            ms_zoom_pic.text(m1_theo_mz + 2.5, max(ms_zoom_bp_i - ms_zoom_offset_i, ms1_pr_i * 0.8),
                             '[M+H2] Isotope score = %.1f' % m2_score,
                             verticalalignment='top', horizontalalignment='right',
                             color='red', fontsize=7)

        # plot the isotope score
        ms_zoom_pic.text(m1_theo_mz + 0.2, ms_zoom_bp_i + 3 * ms_zoom_offset_i,
                         'Isotope score = %.1f' % isotope_score,
                         verticalalignment='top', horizontalalignment='right',
                         color='magenta', fontsize=10)

        # XIC spectrum start
        xic_rt_min = min(xic_rt_lst)
        xic_rt_max = max(xic_rt_lst)
        xic_rt_label_shift = (xic_rt_max - xic_rt_min) * 0.04
        xic_pic.plot(xic_rt_lst, xic_i_lst, alpha=0.7, color='grey')
        xic_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        markerline, stemlines, baseline = xic_pic.stem([ms1_rt], [max(xic_i_lst)], markerfmt=' ')
        plt.setp(stemlines, color='magenta', linewidth=3, alpha=0.3)
        markerline, stemlines, baseline = xic_pic.stem([ms2_rt], [max(xic_i_lst)], '--', markerfmt=' ')
        plt.setp(stemlines, color='blue', linewidth=2, alpha=0.3)
        xic_pic.text(ms1_rt - xic_rt_label_shift, max(xic_i_lst) * 0.98, 'MS', fontsize=8, color='magenta')
        xic_pic.text(ms2_rt, max(xic_i_lst) * 0.98, 'MS/MS', fontsize=8, color='blue')
        xic_pic.set_xlabel("Scan time (min)", fontsize=10, labelpad=-1)
        xic_pic.set_ylabel("Intensity", fontsize=10)
        xic_pic.set_xlim([xic_rt_min, xic_rt_max])

        # prepare DataFrame for msms zoomed plot
        # plot color markers for zoomed MS2 first. Then overlay with zoomed spectra. Plot full ms2 in the last step.

        _msms_max = ms2_df['i'].max()
        _msms_low_delta = _msms_low_df['i'].max() * 0.06
        _msms_high_delta = _msms_high_df['i'].max() * 0.06

        _matched_fa_df = ident_info_dct['MATCHED_FA_INFO']
        _matched_lyso_df = ident_info_dct['MATCHED_LYSO_INFO']
        _other_signals_df = ident_info_dct['OTHER_SIGNALS_INFO']

        overall_score = ident_info_dct['Overall_score']

        ident_col_labels = ('Proposed_structure', 'Score')
        _ident_table_df = pd.DataFrame(data={'Proposed_structure': abbr, 'Score': overall_score}, index=[0])
        ident_table_vals = map(list, _ident_table_df.values)
        ident_col_width_lst = [0.6, 0.15]
        ident_table = msms_pic.table(cellText=ident_table_vals, colWidths=ident_col_width_lst,
                                     colLabels=ident_col_labels, loc='upper center', cellLoc='center')
        ident_table.set_fontsize(8)

        if _matched_fa_df.shape[0] > 0:
            matched_fa_mz_lst = _matched_fa_df['mz'].tolist()
            for _i_fa, _fa_se in _matched_fa_df.iterrows():

                markerline, stemlines, baseline = msms_pic.stem([_fa_se['mz']], [_fa_se['i']], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3)
                if _fa_se['mz'] <= 400:
                    markerline, stemlines, baseline = msms_low_pic.stem([_fa_se['mz']], [_fa_se['i']], markerfmt=' ')
                    plt.setp(stemlines, color='#00ccff', linewidth=3)
                    _msms_low_peak_y = float(_fa_se['i'])
                    _msms_low_peak_str = '%.4f' % _fa_se['mz']
                    msms_low_pic.text(_fa_se['mz'], _msms_low_peak_y, _msms_low_peak_str, fontsize=6, color='#00ccff')
                    msms_low_pic.text(_fa_se['mz'], _msms_low_peak_y + _msms_low_delta, _fa_se['Proposed_structures'],
                                      color='#00ccff', fontsize=8)
                else:
                    markerline, stemlines, baseline = msms_high_pic.stem([_fa_se['mz']], [_fa_se['i']], markerfmt=' ')
                    plt.setp(stemlines, color='#00ccff', linewidth=3)
                    _msms_high_peak_y = float(_fa_se['i'])
                    _msms_high_peak_str = '%.4f' % _fa_se['mz']
                    msms_high_pic.text(_fa_se['mz'], _msms_high_peak_y, _msms_high_peak_str, fontsize=6, color='#00ccff')
                    msms_high_pic.text(_fa_se['mz'], _msms_high_peak_y + _msms_high_delta, _fa_se['Proposed_structures'],
                                       color='#00ccff', fontsize=8)
        else:
            matched_fa_mz_lst = []

        if _matched_lyso_df.shape[0] > 0:
            matched_lyso_mz_lst = _matched_lyso_df['mz'].tolist()
            for _i_lyso, _lyso_se in _matched_lyso_df.iterrows():
                markerline, stemlines, baseline = msms_pic.stem([_lyso_se['mz']], [_lyso_se['i']], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3)
                markerline, stemlines, baseline = msms_high_pic.stem([_lyso_se['mz']], [_lyso_se['i']], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3)
                _msms_high_peak_y = float(_lyso_se['i'])
                _msms_high_peak_str = '%.4f' % _lyso_se['mz']
                msms_high_pic.text(_lyso_se['mz'], _msms_high_peak_y, _msms_high_peak_str, fontsize=6, color='#00ccff')
                msms_high_pic.text(_lyso_se['mz'], _msms_high_peak_y + _msms_high_delta, _lyso_se['Proposed_structures'],
                                   color='#00ccff', fontsize=8)
        else:
            matched_lyso_mz_lst = []

        if _other_signals_df.shape[0] > 0:
            for _i_other_sig, _other_sig_se in _other_signals_df.iterrows():
                if _other_sig_se['mz'] in matched_lyso_mz_lst or _other_sig_se['mz'] in matched_fa_mz_lst:
                    pass
                else:
                    markerline, stemlines, baseline = msms_pic.stem([_other_sig_se['mz']], [_other_sig_se['i']],
                                                                    markerfmt=' ')
                    plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                    if _other_sig_se['mz'] <= 400:
                        markerline, stemlines, baseline = msms_low_pic.stem([_other_sig_se['mz']], [_other_sig_se['i']],
                                                                            markerfmt=' ')
                        plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                        _msms_low_peak_str = '%.4f' % _other_sig_se['mz']
                        _msms_low_peak_y = float(_other_sig_se['i'])
                        msms_low_pic.text(_other_sig_se['mz'], _msms_low_peak_y, _msms_low_peak_str, fontsize=6,
                                          color='red')
                    else:
                        markerline, stemlines, baseline = msms_high_pic.stem([_other_sig_se['mz']], [_other_sig_se['i']],
                                                                             markerfmt=' ')
                        plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                        _msms_high_peak_str = '%.4f' % _other_sig_se['mz']
                        _msms_high_peak_y = float(_other_sig_se['i'])
                        msms_high_pic.text(_other_sig_se['mz'], _msms_high_peak_y, _msms_high_peak_str,
                                           fontsize=6, color='red')

        # msms spectrum zoomed < 400 start
        if _msms_low_df.shape[0] > 0:
            msms_low_pic.stem(_msms_low_df['mz'].tolist(),
                              _msms_low_df['i'].tolist(),
                              'black', lw=2, markerfmt=' ')
            msms_low_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_low_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            msms_low_pic.set_ylabel("Intensity", fontsize=10)
            msms_low_pic.set_xlim([min(_msms_low_df['mz'].tolist()) - 1, 400])
            msms_low_pic.set_ylim([0, max(_msms_low_df['i'].tolist()) * 1.3])
        else:
            pass

        # msms spectrum zoomed > 400 start
        if _msms_high_df.shape[0] > 0:
            msms_high_pic.stem(_msms_high_df['mz'].tolist(),
                               _msms_high_df['i'].tolist(),
                               'black', lw=4, markerfmt=' ')
            msms_high_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            msms_high_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
            msms_high_pic.set_ylabel("Intensity", fontsize=10)
            msms_high_pic.set_xlim([400, ms2_pr_mz + 5])
            msms_high_pic.set_ylim([0, max(_msms_high_df['i'].tolist()) * 1.3])

            # add annotations
            _top_msms_high_df = _msms_high_df.sort_values(by='i', ascending=False)
            _top_msms_high_df = _top_msms_high_df.head(10)
            _msms_high_peak_list = zip(_top_msms_high_df['mz'].tolist(),
                                       _top_msms_high_df['i'].tolist())
            for _msms_high_peak in _msms_high_peak_list:
                _msms_high_peak_str = '%.4f' % _msms_high_peak[0]
                _msms_high_peak_y = _msms_high_peak[1]
                msms_high_pic.text(_msms_high_peak[0], _msms_high_peak_y, _msms_high_peak_str, fontsize=6)
        else:
            pass

        # add specific ion info

        if 'OTHER_FRAG' in specific_check_dct.keys():
            other_frag_df = specific_check_dct['OTHER_FRAG']
            for _idx, _frag_se in other_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i = _frag_se['i']
                _frag_class = _frag_se['LABEL']
                _frag_i_x = min([_frag_i * 5, 0.3 * max(_msms_low_df['i'].tolist())])
                _frag_i_r = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
                markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i_r], markerfmt=' ')
                plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                markerline, stemlines, baseline = msms_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
                plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, fontsize=8, color='red')

        if 'OTHER_NL' in specific_check_dct.keys():
            other_nl_df = specific_check_dct['OTHER_NL']
            for _idx, _nl_se in other_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                _nl_i_x = min([_nl_i * 5, 0.3 * max(_msms_high_df['i'].tolist())])
                _nl_i_r = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
                markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i_r], markerfmt=' ')
                plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                markerline, stemlines, baseline = msms_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
                plt.setp(stemlines, color='red', linewidth=3, alpha=0.5)
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, fontsize=8, color='red')

        if 'TARGET_FRAG' in specific_check_dct.keys():
            target_frag_df = specific_check_dct['TARGET_FRAG']
            for _idx, _frag_se in target_frag_df.iterrows():
                _frag_mz = _frag_se['mz']
                _frag_i = _frag_se['i']
                _frag_class = _frag_se['LABEL']
                _frag_i_x = min([_frag_i * 5, 0.5 * max(_msms_low_df['i'].tolist())])
                _frag_i_r = sorted([max(_msms_low_df['i'].tolist()) * 1.1, _frag_i * 1.1, _frag_i_x])[1]
                markerline, stemlines, baseline = msms_low_pic.stem([_frag_mz], [_frag_i_r], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3, alpha=0.5)
                markerline, stemlines, baseline = msms_pic.stem([_frag_mz], [_frag_i], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3, alpha=0.5)
                msms_low_pic.text(_frag_mz, _frag_i_r, _frag_class, fontsize=8, color='#00ccff')

        if 'TARGET_NL' in specific_check_dct.keys():
            target_nl_df = specific_check_dct['TARGET_NL']
            for _idx, _nl_se in target_nl_df.iterrows():
                _nl_mz = _nl_se['mz']
                _nl_i = _nl_se['i']
                _nl_class = _nl_se['LABEL']
                _nl_i_x = min([_nl_i * 5, 0.3 * max(_msms_low_df['i'].tolist())])
                _nl_i_r = sorted([max(_msms_high_df['i'].tolist()) * 1.1, _nl_i * 1.1, _nl_i_x])[1]
                markerline, stemlines, baseline = msms_high_pic.stem([_nl_mz], [_nl_i_r], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3, alpha=0.5)
                markerline, stemlines, baseline = msms_pic.stem([_nl_mz], [_nl_i], markerfmt=' ')
                plt.setp(stemlines, color='#00ccff', linewidth=3, alpha=0.5)
                msms_high_pic.text(_nl_mz, _nl_i_r, _nl_class, fontsize=8, color='#00ccff')

        # msms spectrum start
        if plot_msp == 1:
            fp_max_i = min(msp_info['rev_abs_i'].tolist())
        else:
            fp_max_i = _msms_max * -0.5
        if len(obs_fp) > 0:
            markerline, stemlines, baseline = msms_pic.stem(obs_fp, [fp_max_i * 1.1] * len(obs_fp),
                                                            ':', markerfmt=' ')
            plt.setp(stemlines, color='#88ff88', alpha=0.4, lw=1)
        if len(missed_fp) > 0:
            markerline, stemlines, baseline = msms_pic.stem(missed_fp, [fp_max_i] * len(missed_fp),
                                                            ':', markerfmt=' ')
            plt.setp(stemlines, color='#999999', alpha=0.4, lw=1)

        min_msp_i = fp_max_i
        min_msp_i_fp = min_msp_i * 1.55
        min_msp_i_fp_h = abs(min_msp_i * 0.5)

        # plot missed ones, so they can be overlaid by the identified ones
        for missed_fp_mz in missed_fp:
            missed_fp_box = patches.Rectangle((missed_fp_mz - 1.75, min_msp_i_fp), 3.5, min_msp_i_fp_h,
                                              facecolor='#999999', edgecolor='none')
            msms_pic.add_patch(missed_fp_box)

        for obs_fp_mz in obs_fp:
            obs_fp_box = patches.Rectangle((obs_fp_mz - 1.75, min_msp_i_fp), 3.5, min_msp_i_fp_h,
                                           facecolor='#88ff88', edgecolor='none')
            msms_pic.add_patch(obs_fp_box)

        if plot_msp == 1:

            msms_pic.stem(msp_info['mz'].tolist(), msp_info['rev_abs_i'].tolist(), '#ff6600',
                          markerfmt=' ', basefmt='k-', zorder=2)
        msms_pic.stem(ms2_df['mz'].tolist(), ms2_df['i'].tolist(), 'black', markerfmt=' ', basefmt='k-', zorder=10)
        msms_pic.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        msms_pic.set_xlabel("m/z", fontsize=10, labelpad=-1)
        msms_pic.set_ylabel("Intensity", fontsize=10)
        msms_pic.set_xlim([min(ms2_df['mz'].tolist()) - 1, ms2_pr_mz + 5])
        msms_pic.set_ylim([min_msp_i * 1.55, _msms_max * 1.6])

        # set title
        xic_title_str = 'XIC of m/z %.4f | @ m/z %.4f ppm=%.2f' % (ms1_pr_mz, lib_mz, ms1_pr_ppm)
        ms_title_str = 'MS @ %.3f min' % ms1_rt
        ms_zoom_title_str = 'Theoretical isotopic distribution for %s %s' % (formula_charged, charge)
        msms_title_str = ('MS/MS for m/z %.4f | DDA rank %d @ %.3f min' % (ms2_pr_mz, func_id, ms2_rt))
        msms_low_str = 'MS/MS zoomed below m/z 400'
        msms_high_str = 'MS/MS zoomed above m/z 400'

        xic_pic.set_title(xic_title_str, color='b', fontsize=10, y=0.98)
        ms_pic.set_title(ms_title_str, color='b', fontsize=10, y=0.98)
        ms_zoom_pic.set_title(ms_zoom_title_str, color='b', fontsize=10, y=0.98)
        msms_pic.set_title(msms_title_str, color='b', fontsize=10, y=0.98)
        msms_low_pic.set_title(msms_low_str, color='b', fontsize=10, y=0.98)
        msms_high_pic.set_title(msms_high_str, color='b', fontsize=10, y=0.98)

        print ('>>> >>> >>> try to plot >>> >>> >>>')
        try:
            plt.savefig(save_img_as, type=img_type, dpi=dpi)
            print ('=====> Image saved as: %s' % save_img_as)
            img_name_end = ''
        except IOError:
            plt.savefig('%s-2.%s' % (save_img_as[:-4], img_type), type=img_type, dpi=dpi)
            print ('=====> Image saved as: %s' % save_img_as)
            img_name_end = '-2'
        plt.close()
        isotope_checker = 0
        return isotope_checker, isotope_score, img_name_end

    else:
        img_name_end = ''
        isotope_checker = 1
        return isotope_checker, isotope_score, img_name_end
