# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import natsort
from FAtoLPP import FAtoLPP
from lpplibs.ExactMassCalc import Elem2Mass
from PLclassDict import PL_Class_dct

# s_hg_txt = '[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])('
# hg_a_s = '[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])(OC(C)=O)COC(C)=O)=O'
# hg_mol = Chem.MolFromSmiles(hg_a_s)
# AllChem.Compute2DCoords(hg_mol)
#
# s1 = r'OC(CCCCCCCCCCCCCCC)=O'
# s2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
# # s2 = r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'

# PC hg
s_hg_txt = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])('
hg_a_s = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(C)=O)COC(C)=O)=O'

# PE hg
# s_hg_txt = 'OP(OCCN)(OC[C@]([H])('
# hg_a_s = 'OP(OCCN)(OC[C@]([H])(OC(C)=O)COC(C)=O)=O'

hg_mol = Chem.MolFromSmiles(hg_a_s)
AllChem.Compute2DCoords(hg_mol)

s1 = r'OC(CCCCCCCCCCCCCCC)=O'
# s2 = r'OC(CCCCCCC/C=C\CCCCCCCC)=O'
# s2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
# s2 = r'OC(CCCCCCC/C=C\C/C=C\C/C=C\CC)=O'
s2 = r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
# s2 = r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CC)=O'
# s2 = r'OC(CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CC)=O'

pl_name = 'oxPAPC_OH_KETO_sameDB2'
pl_name_sdf = pl_name + '.sdf'
pl_name_csv = pl_name + '.csv'
pl_name_png = pl_name + '.png'

sdf2img = 0

s1_obj = FAtoLPP(s1)
s1_lst = s1_obj.get_lpp_all()
s2_obj = FAtoLPP(s2)
s2_lst = s2_obj.get_lpp_all()

# print s1_lst
# print s2_lst

pl_smiles_lst = []
pl_d_lst = []

for _s1 in s1_lst:
    for _s2 in s2_lst:
        _s_s1 = _s1[2]
        _s_s2 = _s2[2]
        tmp_s_pl_lst = [s_hg_txt, _s_s2, ')C', _s_s1, ')=O']
        _s_pl = ''.join(tmp_s_pl_lst)
        tmp_d_pl_lst = ['oxPC(', _s1[0], '/', _s2[0], ')']
        _d_pl = ''.join(tmp_d_pl_lst)

        pl_smiles_lst.append(_s_pl)
        pl_d_lst.append(_d_pl)

sum_lst = zip(pl_d_lst, pl_smiles_lst)

pl_mol_lst = []
for _f in sum_lst:
    tmp_z_mol = Chem.MolFromSmiles(_f[1])
    AllChem.Compute2DCoords(tmp_z_mol)
    AllChem.GenerateDepictionMatching2DStructure(tmp_z_mol, hg_mol)
    tmp_z_mol.SetProp("_Name", _f[0])
    tmp_z_mol.SetProp("LM_ID", _f[0])
    tmp_z_mol.SetProp("COMMON_NAME", _f[0])
    tmp_z_mol.SetProp("SMILES", _f[1])
    pl_mol_lst.append(tmp_z_mol)

if sdf2img == 1:
    img = Draw.MolsToGridImage(pl_mol_lst, molsPerRow=8, legends=pl_d_lst, subImgSize=(400, 400))
    img.save(pl_name_png)
else:
    pass


w = Chem.SDWriter(pl_name_sdf)
mzcalc = Elem2Mass()
exact_mass_lst = []
elem_lst = []
mz_H_lst = []
mz_Na_lst = []

for _mol in pl_mol_lst:
    _smiles = _mol.GetProp("SMILES")
    elem_db = {}
    smiles_lst = list(_smiles)
    elem_db['C'] = smiles_lst.count('C')
    elem_db['O'] = smiles_lst.count('O')
    elem_db['P'] = smiles_lst.count('P')
    elem_db['N'] = smiles_lst.count('N')
    elem_db['dbe'] = smiles_lst.count('=')
    # elem_db['D'] = 9
    elem_db['H'] = smiles_lst.count('C') * 2 + 2 + 4 - 2 * smiles_lst.count('=')  # 9D

    GP_elem_str = ''
    GP_elem_idx_lst = ['C', 'H', 'O', 'P', 'N']  # , 'D'

    for tmp_elem in GP_elem_idx_lst:
        if tmp_elem in elem_db.keys():
            tmp_elem_num = elem_db[tmp_elem]
            tmp_info = tmp_elem + str(tmp_elem_num)
            GP_elem_str += tmp_info

        else:
            pass

    formula_org = mzcalc.get_elem(GP_elem_str)
    exact_mass = mzcalc.get_mass(formula_org)
    formula_H = mzcalc.get_elem(GP_elem_str + 'H')
    mz_H = mzcalc.get_mass(formula_H)
    formula_Na = mzcalc.get_elem(GP_elem_str + 'Na')
    mz_Na = mzcalc.get_mass(formula_Na)
    exact_mass_lst.append(exact_mass)
    mz_H_lst.append(mz_H)
    mz_Na_lst.append(mz_Na)
    elem_lst.append(GP_elem_str)

    # _mol.SetProp('Formula', GP_elem_str)
    # _mol.SetProp('Exact_Mass', str(exact_mass))
    w.write(_mol)
w.close()

info_dct = {'PL_Abbr': pl_d_lst, 'Elem': elem_lst, 'SMILES': pl_smiles_lst, 'ExactMass': exact_mass_lst,
            '[M+H]+': mz_H_lst, '[M+Na]+': mz_Na_lst}

df = pd.DataFrame(data=info_dct, columns=['PL_Abbr', 'Elem', 'ExactMass', '[M+H]+', '[M+Na]+', 'SMILES'])
print df.head()

df.to_csv(pl_name_csv)

print 'fin!'