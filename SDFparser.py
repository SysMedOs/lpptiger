# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from rdkit import Chem
from ExactMassCalc import Elem2Mass

w = Chem.SDWriter('merged_info3.sdf')
suppl = Chem.SDMolSupplier('merged_test.sdf')
mzcalc = Elem2Mass()

# exact_mass_lst = []
# elem_lst = []
# mz_H_lst = []
# mz_Na_lst = []

for m in suppl:
    smiles = Chem.MolToSmiles(m)

    m.SetProp('SMILES', smiles)
    print smiles

    elem_db = {}
    smiles_lst = list(smiles)
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

    # formula_H = mzcalc.get_elem(GP_elem_str + 'H')
    # mz_H = mzcalc.get_mass(formula_H)
    # formula_Na = mzcalc.get_elem(GP_elem_str + 'Na')
    # mz_Na = mzcalc.get_mass(formula_Na)
    # exact_mass_lst.append(exact_mass)
    # mz_H_lst.append(mz_H)
    # mz_Na_lst.append(mz_Na)
    # elem_lst.append(GP_elem_str)

    m.SetProp('Formula', GP_elem_str)
    m.SetProp('Exact_Mass', str(exact_mass))

    w.write(m)

w.close()

print 'SDF created! '

# info_dct = {'PL_Abbr': d_z_lst, 'Elem': elem_lst, 'SMILES': s_z_lst, 'ExactMass': exact_mass_lst,
#                 '[M+H]+': mz_H_lst, '[M+Na]+': mz_Na_lst}
#
# df = pd.DataFrame(data=info_dct, columns=['PL_Abbr', 'Elem', 'ExactMass', '[M+H]+', '[M+Na]+', 'SMILES'])
# print df.head()
#
# df.to_csv(csv_name)

print 'Finished !'
