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
import platform
import ConfigParser as configparser
from numba import int32, float32, float64, vectorize


# temp fix for numba vectorize issue of target = 'parallel' for Linux
if platform.system() == 'Windows':
    print('Parallel processing mode: Windows')
    para_target = 'parallel'
elif platform.system() == 'Linux':
    print('Parallel processing mode: Linux')
    para_target = 'cpu'
else:
    print('Parallel processing mode: other systems')
    para_target = 'cpu'

# setup weight factor
# load configurations
config = configparser.ConfigParser()
config.read('configure.ini')
if config.has_section('settings'):
    user_cfg = 'settings'
else:
    if config.has_section('default'):
        user_cfg = 'default'
    else:
        user_cfg = ''
if len(user_cfg) > 2:
    options = config.options(user_cfg)
    if 'general_mod_lst' in options:
        intensity_factor = float(config.get(user_cfg, 'specsim_m'))
    else:
        intensity_factor = 0.6
    if 'general_mod_lst' in options:
        mz_factor = float(config.get(user_cfg, 'specsim_n'))
    else:
        mz_factor = 0.6
else:
    intensity_factor = 3
    mz_factor = 0.6

print('intensity_factor', intensity_factor, type(intensity_factor), 'mz_factor', mz_factor, type(mz_factor))

@vectorize(([float64(float64, float32)]), target=para_target)
def pr_window_calc_para(mz, delta):
    return mz + delta


@vectorize(([float64(float64, int32)]), target=para_target)
def ppm_window_para(mz, ppm):
    return mz * (1 + 0.000001 * ppm)


@vectorize(([float64(float64, float64)]), target=para_target)
def ppm_calc_para(mz_obs, mz_lib):
    return 1e6 * (mz_obs - mz_lib) / mz_lib


@vectorize(([float64(float64, float64)]), target=para_target)
def wfactor_calc_para(mz, i):
    return (mz ** mz_factor) * (i ** intensity_factor)
