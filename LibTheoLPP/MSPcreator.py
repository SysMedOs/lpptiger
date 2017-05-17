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

# import re
import json


def to_msp(output_obj, lpp_info_dct):

    _lpp_name = ''.join(['Name: ', lpp_info_dct['LM_ID'], ' ', lpp_info_dct['LM_ID']])
    _lpp_id = ''.join(['LM_ID: ', lpp_info_dct['LM_ID']])
    _lpp_pr_dct = json.loads(lpp_info_dct['PRECURSOR_JSON'])
    # '{"[M-H]-": ["%s", %f]}' % (_lpp_formula, _lpp_neg_precursor_mz[0])
    _lpp_pr_typ = ''.join(['Precursor_type: ', _lpp_pr_dct.keys()[0]])
    _lpp_pr_mz = ''.join(['PrecursorMZ: ', str(_lpp_pr_dct[_lpp_pr_dct.keys()[0]][1])])
    _lpp_formula = ''.join(['Formula: ', _lpp_pr_dct[_lpp_pr_dct.keys()[0]][0]])
    _lpp_peaks_dct = json.loads(lpp_info_dct['MSP_JSON'])
    _lpp_num_peaks = ''.join(['Num Peaks: ', str(len(_lpp_peaks_dct.keys()))])

    ion_lst = [_lpp_name, _lpp_id, _lpp_pr_typ, _lpp_pr_mz, _lpp_formula, _lpp_num_peaks]
    for _ion in _lpp_peaks_dct.keys():
        _ion_name = ''.join(['"', _ion, ' | ',  _lpp_peaks_dct[_ion]['formula'], '"'])
        _ion_info_txt = ' '.join([str(_lpp_peaks_dct[_ion]['mz']), str(_lpp_peaks_dct[_ion]['i']), _ion_name])
        ion_lst.append(_ion_info_txt)
    ion_lst.append('\n')

    output_obj.writelines('\n'.join(ion_lst))
