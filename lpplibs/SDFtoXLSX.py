# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function

import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit import RDConfig


def sdf2xlsx(usr_sdf, usr_output):

    usr_df =PandasTools.LoadSDF(usr_sdf)
    print(usr_df.shape)

    usr_export = usr_df.to_excel(usr_output)

    return usr_df

# # not functioning
# f = r'new_method_sdf_PE_max_1keto_1lessDB.sdf'
#
# out = 'new_method_sdf_PE_max_1keto_1lessDB.xlsx'
#
# a = sdf2xlsx(f, out)
#
# print(a.shape)
