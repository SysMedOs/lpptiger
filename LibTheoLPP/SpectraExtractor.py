# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of TheoLPP.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from __future__ import print_function
from __future__ import division

import re

from matplotlib import pyplot as plt
import pandas as pd
import pymzml


class SpectraExtractor(object):

    def __init__(self, mzml_path, mz_range=(350, 1000), rt_range=(0, -1), rt_unit='min',
                 ms1_threshold=1000, ms2_threshold=10,
                 ms1_precision=(100, 'ppm'), ms2_precision=(0.5, 'Da')):
        """
        hot patch for waters mzML fields
        ('MS:1000016', ['value']) --> Start time of the scan
        ('MS:1000744', ['value']) --> precursor for MS2
        ('MS:1000042', ['value']) --> precursor intensity
        ('MS:1000796', ['value']) --> spectrum title
        070120_CM_neg_70min_BOTH_I.1.1. File:"070120_CM_neg_70min_BOTH_I.raw", NativeID:"function=1 process=0 scan=71"
        ('MS:1000514', ['name']) --> m/z array
        ('MS:1000515', ['name']) --> intensity array
        """
        waters_obo_lst = (('MS:1000016', ['value']), ('MS:1000744', ['value']), ('MS:1000042', ['value']),
                          ('MS:1000796', ['value']), ('MS:1000514', ['name']), ('MS:1000515', ['name']))
        for _obo in waters_obo_lst:
            if _obo not in pymzml.minimum.MIN_REQ:
                pymzml.minimum.MIN_REQ.append(_obo)
            else:
                pass
        # end hot patch

        # mz = -1 --> max mz set as 2000 m/z
        if len(mz_range) == 2:
            mz_min = mz_range[0]
            if mz_range[1] == -1:
                mz_max = 2000
            else:
                mz_max = mz_range[1]
        else:
            mz_min = 0
            mz_max = 2000

        # rt = -1 --> max rt set as 1000000 min
        if rt_unit == 's' and len(rt_range) == 2:
            rt_start_min = rt_range[0] / 60
            if rt_range[1] == -1:
                rt_end_min = 1000000
            else:
                rt_end_min = rt_range[1] / 60
        elif rt_unit == 'min' and len(rt_range) == 2:
            rt_start_min = rt_range[0]
            if rt_range[1] == -1:
                rt_end_min = 1000000
            else:
                rt_end_min = rt_range[1]
        else:
            rt_start_min = 0
            rt_end_min = 1000000

        self.mzml = mzml_path
        self.mz_min = mz_min
        self.mz_max = mz_max
        self.rt_start_min = rt_start_min
        self.rt_end_min = rt_end_min
        self.ms1_threshold = ms1_threshold
        self.ms2_threshold = ms2_threshold
        self.ms1_precision = ms1_precision
        self.ms2_precision = ms2_precision

    def get_ms_all(self):

        _usr_spectra = pymzml.run.Reader(self.mzml)
        _pre_ms_df = pd.DataFrame(columns=['mz', 'i', 'rt', 'scan_id'])

        _spec_title_obo = 'MS:1000796'

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)')

        # # set default settings#
        for _spectrum in _usr_spectra:

            if 'MS:1000016' in _spectrum.keys() and self.rt_start_min <= _spectrum['MS:1000016'] <= self.rt_end_min:
                _spectrum_title = _spectrum[_spec_title_obo]
                _function_checker = _function_re.match(_spectrum_title)
                if _function_checker:
                    _function = _function_checker.groups()[2]

                    if _function == '1':
                        print ('Function: {func}, Scan_num: {scan}, Scan_time: {rt};'
                               .format(func=_function, scan=_spectrum['id'], rt=_spectrum['MS:1000016']))
                        _top_peaks_lst = _spectrum.peaks
                        _top_peaks_df = pd.DataFrame(data=_top_peaks_lst, columns=['mz', 'i'])
                        _top_peaks_df['rt'] = _spectrum['MS:1000016']
                        _top_peaks_df['scan_id'] = _spectrum['id']
                        _pre_ms_df = _pre_ms_df.append(_top_peaks_df)
                    else:
                        pass
                else:
                    pass
            else:
                pass
        _ms_df = _pre_ms_df[_pre_ms_df['i'] >= self.ms1_threshold]
        _ms_df = _ms_df.sort_values(by='mz')
        print('_ms_df_shape', _ms_df.shape)
        return _ms_df

    def get_scans(self, dda_top=12):

        _usr_spectra = pymzml.run.Reader(self.mzml)
        _pre_ms_df = pd.DataFrame(columns=['function', 'scan_id', 'rt', 'mz', 'i'])

        _ms2_function_lst = [str(i) for i in range(2, dda_top + 2)]

        _spec_title_obo = 'MS:1000796'

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')

        for _spectrum in _usr_spectra:
            if 'MS:1000016' in _spectrum.keys() and self.rt_start_min <= _spectrum['MS:1000016'] <= self.rt_end_min:
                if _spec_title_obo in _spectrum.keys() and _function_re.match(_spectrum[_spec_title_obo]):
                    _function_checker = _function_re.match(_spectrum[_spec_title_obo])
                    _function = _function_checker.groups()[2]
                    _scan = _function_checker.groups()[5]

                    if _function == '1':
                        print('Function: {func}, Scan_num: {scan}, Scan_time: {rt};'
                              .format(func=_function, scan=_spectrum['id'], rt=_spectrum['MS:1000016']))
                        _top_peaks_lst = _spectrum.peaks
                        _top_peaks_df = pd.DataFrame(data=_top_peaks_lst, columns=['mz', 'i'])
                        _top_peaks_df['function'] = _function
                        _top_peaks_df['rt'] = _spectrum['MS:1000016']
                        _top_peaks_df['scan_id'] = _scan
                        _top_peaks_df = _top_peaks_df[(_top_peaks_df['mz'] >= self.mz_min) &
                                                      (_top_peaks_df['mz'] <= self.mz_max)]
                        _top_peaks_df = _top_peaks_df[_top_peaks_df['i'] >= self.ms1_threshold]
                        _pre_ms_df = _pre_ms_df.append(_top_peaks_df)

                    if _function in _ms2_function_lst:
                        print('Function: {func}, Scan_num: {scan}, Scan_time: {rt};'
                              .format(func=_function, scan=_spectrum['id'], rt=_spectrum['MS:1000016']))
                        _top_peaks_lst = [(_spectrum['MS:1000744'], 0)]  # Get PR info only
                        _top_peaks_df = pd.DataFrame(data=_top_peaks_lst, columns=['mz', 'i'])  # PR i default -> 0.0
                        _top_peaks_df['function'] = _function
                        _top_peaks_df['rt'] = _spectrum['MS:1000016']
                        _top_peaks_df['scan_id'] = _scan
                        _pre_ms_df = _pre_ms_df.append(_top_peaks_df)

                else:
                    print('NOT MS level')
        _scans_df = _pre_ms_df.sort_values(by=['rt'])
        print('_ms_df_shape', _scans_df.shape)
        return _scans_df

    def get_xic(self, mz):
        _usr_spectra = pymzml.run.Reader(self.mzml)
        if self.ms1_precision[1] == 'ppm':
            _mz_min = mz * (1 - self.ms1_precision[0] / 1000000)
            _mz_max = mz * (1 + self.ms1_precision[0] / 1000000)
        elif self.ms1_precision[1].upper() == 'DA' or self.ms1_precision[1].upper() == 'MZ':
            _mz_min = mz - self.ms1_precision[0]
            _mz_max = mz + self.ms1_precision[0]
        else:
            # set default to 100 ppm
            _mz_min = mz * 0.9999
            _mz_max = mz * 1.0001

        _xic_df = pd.DataFrame(columns=['function', 'scan_id', 'rt', 'mz', 'i'])

        _query_code = '{mz_min} <= mz <= {mz_max} & i >= {ms1_th}'.format(mz_min=_mz_min, mz_max=_mz_max,
                                                                          ms1_th=self.ms1_threshold)

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')
        _spec_title_obo = 'MS:1000796'

        for _spectrum in _usr_spectra:
            if 'MS:1000016' in _spectrum.keys() and self.rt_start_min <= _spectrum['MS:1000016'] <= self.rt_end_min:
                if _spec_title_obo in _spectrum.keys() and _function_re.match(_spectrum[_spec_title_obo]):
                    _function_checker = _function_re.match(_spectrum[_spec_title_obo])
                    _function = _function_checker.groups()[2]
                    _scan = _function_checker.groups()[5]

                    if _function == '1':
                        _top_peaks_lst = _spectrum.peaks
                        _top_peaks_df = pd.DataFrame(data=_top_peaks_lst, columns=['mz', 'i'])
                        _top_peaks_df['function'] = _function
                        _top_peaks_df['rt'] = _spectrum['MS:1000016']
                        _top_peaks_df['scan_id'] = _scan
                        _top_peaks_df = _top_peaks_df.query(_query_code)
                        _top_peaks_df = _top_peaks_df.sort_values(by='i', ascending=False)
                        _xic_df = _xic_df.append(_top_peaks_df.head(1))
                        # _xic_df = _xic_df.append(_top_peaks_df)
                    else:
                        pass
                else:
                    pass
            else:
                pass

        _xic_df = _xic_df.sort_values(by='rt')
        print('_ms_df_shape', _xic_df.shape)
        return _xic_df

    def get_ms2(self, mz, dda_top=12):
        _usr_spectra = pymzml.run.Reader(self.mzml)
        if self.ms1_precision[1] == 'ppm':
            _mz_min = mz * (1 - self.ms1_precision[0] / 1000000)
            _mz_max = mz * (1 + self.ms1_precision[0] / 1000000)
        elif self.ms1_precision[1].upper() == 'DA' or self.ms1_precision[1].upper() == 'MZ':
            _mz_min = mz - self.ms1_precision[0]
            _mz_max = mz + self.ms1_precision[0]
        else:
            # set default to 0.5 Da
            _mz_min = mz - 0.5
            _mz_max = mz + 0.5
        _ms2_function_lst = [str(i) for i in range(2, dda_top + 2)]

        _ms2_df = pd.DataFrame(columns=['function', 'scan_id', 'rt', 'mz', 'i'])
        _ms2_rt_lst = []

        _query_code = '{mz_min} <= mz <= {mz_max} & i >= {ms2_th}'.format(mz_min=_mz_min, mz_max=_mz_max,
                                                                          ms2_th=self.ms2_threshold)

        _function_re = re.compile(r'(.*)(function=)(\d{1,2})(.*)(scan=)(\d*)(.*)')
        _spec_title_obo = 'MS:1000796'

        for _spectrum in _usr_spectra:
            if 'MS:1000016' in _spectrum.keys() and self.rt_start_min <= _spectrum['MS:1000016'] <= self.rt_end_min:
                if _spec_title_obo in _spectrum.keys() and _function_re.match(_spectrum[_spec_title_obo]):
                    _function_checker = _function_re.match(_spectrum[_spec_title_obo])
                    _function = _function_checker.groups()[2]
                    _scan = _function_checker.groups()[5]

                    if _function in _ms2_function_lst:
                        print('Function: {func}, Scan_num: {scan}, Scan_time: {rt};'
                              .format(func=_function, scan=_spectrum['id'], rt=_spectrum['MS:1000016']))

                        if _mz_min <= _spectrum['MS:1000744'] <= _mz_max:
                            _top_peaks_lst = _spectrum.peaks
                            _top_peaks_df = pd.DataFrame(data=_top_peaks_lst, columns=['mz', 'i'])
                            _top_peaks_df['function'] = _function
                            _top_peaks_df['rt'] = _spectrum['MS:1000016']
                            _ms2_rt_lst.append(_spectrum['MS:1000016'])
                            _top_peaks_df['scan_id'] = _scan
                            _ms2_df = _ms2_df.append(_top_peaks_df)
                    else:
                        pass
                else:
                    pass
            else:
                pass

        _ms2_df = _ms2_df.sort_values(by='rt')
        print('_ms_df_shape', _ms2_df.shape)
        return _ms2_df, _ms2_rt_lst

if __name__ == '__main__':

    mzml = r'D:\Synapt_rawspectra\oxPLstd\MostAbs10\180816_oxPE_10ng.mzML'

    spectra_obj = SpectraExtractor(mzml, mz_range=(350, 1000), rt_range=(3, 20),
                                   ms1_threshold=1000, ms2_threshold=10, ms1_precision=(100, 'ppm'),
                                   ms2_precision=(0.5, 'Da'))
    # ms_df = spectra_obj.get_ms_all()
    # ms_df = spectra_obj.get_scans()
    xic_df = spectra_obj.get_xic(716.522400)
    ms2_df, ms2_rt_lst = spectra_obj.get_ms2(716.522400)

    fig, axarr = plt.subplots(len(ms2_rt_lst) + 1)
    axarr[0].plot(xic_df['rt'], xic_df['i'], alpha=0.3)
    for _ms2_rt in ms2_rt_lst:
        _idx = ms2_rt_lst.index(_ms2_rt) + 1
        _tmp_ms2_df = ms2_df[ms2_df['rt'] == _ms2_rt]
        axarr[_idx].stem(_tmp_ms2_df['mz'], _tmp_ms2_df['i'], lw=2, alpha=0.3)
    plt.show()
