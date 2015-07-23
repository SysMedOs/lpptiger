# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import os
from rdkit import Chem

f = os.getcwd() + '/SDFs'
print f
os.chdir(f)
print os.listdir(f)
sdf_lst = [x for x in os.listdir(f) if os.path.isfile(x) and os.path.splitext(x)[1] == '.sdf']
# sdf_lst = ['PC(18-0_18-1(9-Z)).sdf', 'PC(18-0_18-2(9-Z;12-Z)).sdf']
print 'sdf_lst', sdf_lst
w = Chem.SDWriter('merged_test.sdf')

#
# merged_sdf = Chem.SDMolSupplier('merge_test.sdf')
merged_dct = {}
merged_lst = []

for tmp_sdf in sdf_lst:
    print tmp_sdf
    suppl = Chem.SDMolSupplier(tmp_sdf)

    for m in suppl:
        m_smiles = Chem.MolToSmiles(m)
        if m_smiles not in merged_dct.keys():
            merged_dct[m_smiles] = m
            print 'add'
        else:
            print 'Inside already! --> skip!', m_smiles

for _m in merged_dct.keys():
    w.write(merged_dct[_m])
w.close()

print 'finished !'

