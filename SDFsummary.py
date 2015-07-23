# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import pandas as pd
from rdkit import Chem
from ExactMassCalc import Elem2Mass


mzcalc = Elem2Mass()

description_lst = []
smiles_lst = []
formula_lst = []
exact_mass_lst = []
mz_H_lst = []
mz_Na_lst = []
mz_K_lst = []

# Get the SDF file
suppl = Chem.SDMolSupplier('merged_info3.sdf')

print 'Got SDF, start processing ...'

for m in suppl:
    _description = m.GetProp('COMMON_NAME')
    _smiles = m.GetProp('SMILES')
    _formula = m.GetProp('Formula')
    _exact_mass = m.GetProp('Exact_Mass')

    # Get parsed formula
    _formula_H = mzcalc.get_elem(_formula + 'H')
    _formula_Na = mzcalc.get_elem(_formula + 'Na')
    _formula_K = mzcalc.get_elem(_formula + 'K')

    _mz_H = mzcalc.get_mass(_formula_H)
    _mz_Na = mzcalc.get_mass(_formula_Na)
    _mz_K = mzcalc.get_mass(_formula_K)

    exact_mass_lst.append(_exact_mass)
    mz_H_lst.append(_mz_H)
    mz_Na_lst.append(_mz_Na)
    mz_K_lst.append(_mz_K)

    description_lst.append(_description)
    smiles_lst.append(_smiles)
    formula_lst.append(_formula)
    print ' finished calculation of: ', _description

print 'Finished reading SDF! ======> Generating summary info...'

info_dct = {'PL_Abbr': description_lst, 'Elem': formula_lst, 'SMILES': smiles_lst, 'ExactMass': exact_mass_lst,
            '[M+H]+': mz_H_lst, '[M+Na]+': mz_Na_lst, '[M+K]+': mz_K_lst}

df = pd.DataFrame(data=info_dct, columns=['PL_Abbr', 'Elem', 'ExactMass', '[M+H]+', '[M+Na]+', '[M+K]+', 'SMILES'])
print df.head()
print df.tail()

df.to_csv('SDF_summary_table_lite.csv')

print 'Summary CSV created!'