# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
import natsort
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PLparser import FAabbr, GPabbr
from FAox import FAox
from ExactMassCalc import Elem2Mass


class PLox:

    def __init__(self):
        pass

    def gen_ox_smiles(self, pl_abbr, isotopelabel=None):

        # Use re groups to match PL expressions
        # e.g. PC(18:2(9-Z;12-Z)/18:1(9-Z))
        # (u'PC', u'(', u'18:2(9-Z;12-Z)', u'/', u'18:1(11-Z)', u')')
        pl_std = re.compile('(P[ACEGSI])([(])([\w\-;:()]*)(/)([\w\-;:()]*)([)])')
        pl_checker = pl_std.match(pl_abbr)
        if pl_checker:
            pl_lst = pl_checker.groups()
            hg_abbr = pl_lst[0]
            sn1_abbr = pl_lst[2]
            sn2_abbr = pl_lst[4]
            sn1 = FAabbr(sn1_abbr)
            sn2 = FAabbr(sn2_abbr)

            s_sn1_lst = sn1.FAstr2SMILES(exact_pos=1)
            s_sn2_lst = sn2.FAstr2SMILES(exact_pos=1)
            s_sn1_txt = s_sn1_lst[0]
            s_sn2_txt = s_sn2_lst[0]

            if hg_abbr == 'PA':
                s_hg_txt = 'OP(O)(OC[C@]([H])('

            elif hg_abbr == 'PC':
                if isotopelabel.lower() == 'd9':
                    s_hg_txt = '[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])('
                else:
                    s_hg_txt = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])('

            elif hg_abbr == 'PE':
                s_hg_txt = 'OP(OCC[N])(OC[C@]([H])('

            elif hg_abbr == 'PG':
                s_hg_txt = 'OP(OCC(O)CO)(OC[C@]([H])('

            elif hg_abbr == 'PS':
                s_hg_txt = 'OP(OCC(C(O)=O)N)(OC[C@]([H])('
            else:
                # as PA
                s_hg_txt = 'OP(O)(OC[C@]([H])('

            # the order is: HG, sn2, sn1 e.g. '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(sn2)=O)COC(sn1)=O)=O'
            tmp_s_pl_lst = [s_hg_txt, s_sn2_txt, ')C', s_sn1_txt, ')=O']
            s_pl = ''.join(tmp_s_pl_lst)

            s_pl_lst = [s_pl]
            d_pl_lst = [pl_abbr]

            ox_sn1 = FAox(sn1_abbr)
            print 's_sn1_txt', s_sn1_txt
            (s_sn1_lst, d_sn1_lst) = ox_sn1.gen_ox_smiles(s_sn1_txt)
            ox_sn2 = FAox(sn2_abbr)
            (s_sn2_lst, d_sn2_lst) = ox_sn2.gen_ox_smiles(s_sn2_txt)

            if len(s_sn1_lst) == 0:
                s_sn1_lst = [s_sn1_txt]
                d_sn1_lst = [sn1_abbr]
            if len(s_sn2_lst) == 0:
                s_sn1_lst = [s_sn2_txt]
                d_sn1_lst = [sn2_abbr]
            sn1_info_lst = zip(s_sn1_lst, d_sn1_lst)
            sn2_info_lst = zip(s_sn2_lst, d_sn2_lst)
            natsort.natsorted(sn1_info_lst, key=lambda t: t[1])
            natsort.natsorted(sn2_info_lst, key=lambda t: t[1])

            s_pl_lst = []
            d_pl_lst = []

            for tmp_sn1 in sn1_info_lst:
                for tmp_sn2 in sn2_info_lst:
                    # the order is: HG, sn2, sn1 e.g.
                    tmp_s_ox_pl_lst = [s_hg_txt, tmp_sn2[0], ')C', tmp_sn1[0], ')=O']
                    tmp_s_ox_pl = ''.join(tmp_s_ox_pl_lst)
                    s_pl_lst.append(tmp_s_ox_pl)
                    tmp_d_ox_pl = hg_abbr + '(' + tmp_sn1[1] + '/' + tmp_sn2[1] + ')'
                    d_pl_lst.append(tmp_d_ox_pl)

        return s_pl_lst, d_pl_lst

# s_x_txt = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)COC(CCCCCCC/C=C\\CCCCCCCC)=O)=O'
# x = 'PC(16:0/18:2(9-Z;12-Z))'
# x = 'PC(16:0/20:4(5-Z;8-Z;11-Z;14-Z))'
# lst_obj = open('list.csv')
# lst = []
# for a in lst_obj.readlines():
#     a = a.strip('\n')
#     lst.append(a)

lst = ['PC(16:0/18:2(9-Z;12-Z))']
print lst

# hg = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(C)=O)COC(C)=O)=O'
hg = '[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])(OC(C)=O)COC(C)=O)=O'
hg_mol = Chem.MolFromSmiles(hg)
AllChem.Compute2DCoords(hg_mol)

for x in lst:
    # try:
    print x
    i = '[d9]PC(16-0_18-2(9-Z;12-Z))'
    img_name = str(i) + '.png'
    sdf_name = str(i) + '.sdf'
    csv_name = str(i) + '.csv'
    p = PLox()

    (s_z_lst, d_z_lst) = p.gen_ox_smiles(x,isotopelabel='d9')
    # print s_z_lst

    z_mol_lst = []
    z_info_lst = zip(s_z_lst, d_z_lst)
    for tmp_z in z_info_lst:
        tmp_d_z = tmp_z[1]
        # print tmp_z[0], tmp_d_z
        tmp_z_mol = Chem.MolFromSmiles(tmp_z[0])
        AllChem.Compute2DCoords(tmp_z_mol)
        AllChem.GenerateDepictionMatching2DStructure(tmp_z_mol, hg_mol)
        tmp_z_mol.SetProp("_Name", tmp_d_z)
        tmp_z_mol.SetProp("LM_ID", tmp_d_z)
        tmp_z_mol.SetProp("COMMON_NAME", tmp_d_z)
        z_mol_lst.append(tmp_z_mol)
    #
    # row_num = len(z_mol_lst)//2
    # if 1 <= row_num <= 9:
    #     pass
    # else:
    #     row_num = 9

    # img = Draw.MolsToGridImage(z_mol_lst, molsPerRow=7, legends=d_z_lst, subImgSize=(300, 300))
    # img = Draw.MolToImage(z_mol, size=(400, 400))
    # img.save(img_name)
    # img.show()

    w = Chem.SDWriter(sdf_name)
    for m in z_mol_lst:
        w.write(m)
    w.close()

    mzcalc = Elem2Mass()
    exact_mass_lst = []
    elem_lst = []
    mz_H_lst = []
    mz_Na_lst = []
    for _smiles in s_z_lst:
        elem_db = {}
        smiles_lst = list(_smiles)
        elem_db['C'] = smiles_lst.count('C')
        elem_db['O'] = smiles_lst.count('O')
        elem_db['P'] = smiles_lst.count('P')
        elem_db['N'] = smiles_lst.count('N')
        elem_db['dbe'] = smiles_lst.count('=')
        elem_db['H'] = smiles_lst.count('C') * 2 + 2 + 4 - 2 * smiles_lst.count('=')

        GP_elem_str = ''
        GP_elem_idx_lst = ['C', 'H', 'O', 'P', 'N']

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

    info_dct = {'PL_Abbr': d_z_lst, 'Elem': elem_lst, 'SMILES': s_z_lst, 'ExactMass': exact_mass_lst,
                '[M+H]+': mz_H_lst, '[M+Na]+': mz_Na_lst}

    df = pd.DataFrame(data=info_dct, columns=['PL_Abbr', 'Elem', 'ExactMass', '[M+H]+', '[M+Na]+', 'SMILES'])
    print df.head()

    df.to_csv(csv_name)

    print 'SDF created!'
    #
    # except:
    #     pass
