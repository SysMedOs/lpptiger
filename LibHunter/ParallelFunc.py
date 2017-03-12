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

from numba import int32, float32, float64, vectorize


@vectorize(([float64(float64, float32)]), target='parallel')
def pr_window_calc_para(mz, delta):
    return mz + delta


@vectorize(([float64(float64, int32)]), target='parallel')
def ppm_window_para(mz, ppm):
    return mz * (1 + 0.000001 * ppm)


@vectorize(([float64(float64, float64)]), target='parallel')
def ppm_calc_para(mz_obs, mz_lib):
    return 1e6 * (mz_obs - mz_lib) / mz_lib
