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
import ConfigParser as configparser
from numba import int32, float32, float64, vectorize


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

    # numba vectorize issue of target = 'parallel' for PC without CUDA compatible GPU
    if 'parallel_target' in options:
        if config.get(user_cfg, 'parallel_target') == 'CPU':
            para_target = 'cpu'
            print('Parallel processing mode -> CPU')
        elif config.get(user_cfg, 'parallel_target') in ['CPU_and_GPU', 'GPU', 'CPUandGPU', 'CPUGPU', 'parallel']:
            para_target = 'parallel'
            print('Parallel processing mode -> CPU_and_GPU')
        else:
            para_target = 'cpu'
            print('Parallel processing mode -> CPU')
    else:
        para_target = 'cpu'
        print('Parallel processing mode -> CPU')

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
    para_target = 'cpu'

print('Weight factor for Similarity score: intensity_factor', intensity_factor, 'mz_factor', mz_factor)


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
