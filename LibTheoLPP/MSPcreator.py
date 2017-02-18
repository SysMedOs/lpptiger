# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

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
        _ion_name = ''.join(['"', _ion, ' | ',  _lpp_peaks_dct[_ion][2], '"'])
        _ion_info_txt = ' '.join([str(_lpp_peaks_dct[_ion][0]), str(_lpp_peaks_dct[_ion][1]), _ion_name])
        ion_lst.append(_ion_info_txt)
    ion_lst.append('\n')

    output_obj.writelines('\n'.join(ion_lst))
