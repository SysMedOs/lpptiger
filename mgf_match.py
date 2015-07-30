# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import pandas as pd
from pyteomics import mgf

csv = 'd9_oxPLPC.csv'
mgf_in = '150715_d9_oxPLPC_100ng_IIIpos.mgf'
mgf_out = '150715_d9_oxPLPC_100ng_IIIpos_lite.mgf'

ppm = 20

df = pd.read_csv(csv, index_col=0)

print df.head()
mz_lst = []
mz_H_lst = df['[M+H]+_d9'].tolist()
mz_Na_lst = df['[M+Na]+_d9'].tolist()
mz_lst.extend(mz_H_lst)
mz_lst.extend(mz_Na_lst)

mz_window_lst = []
for _mz in mz_lst:
    mz_window_lst.append(((_mz * (1 - 0.000001 * ppm)), (_mz * (1 + 0.000001 * ppm))))

# print mz_window_lst

spectra = mgf.read(mgf_in)

reduced_spectra = []

for _spec in spectra:
    _header_dct = _spec['params']
    _pep_info_lst = _header_dct['pepmass']
    _pep_mz = float(_pep_info_lst[0])
    _pep_i = _pep_info_lst[1]
    # print _pep_mz
    try:
        _pep_chg_lst = _header_dct['charge']
        _pep_chg = _pep_chg_lst[0]
    except KeyError:
        _pep_chg = 'X'

    spec_checker = 0

    for _mz_ref in mz_window_lst:

        # if _mz_ref[0] < _pep_mz < _mz_ref[1]:
        if 766 < _pep_mz < 768:
            spec_checker += 1
            _mz_ref_x = _mz_ref
        else:
            pass
    if spec_checker > 0:
        reduced_spectra.append(_spec)
        print 'Fit====================> ', _pep_mz
    else:
        # print 'Not fit, pass...'
        pass
w_spec = mgf.write(reduced_spectra, output=mgf_out)

print 'fin!'



