# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import time
import json

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from lpplibs.PLParser import PLParser
from lpplibs.DBoxTheo import fa_link_filter, oxidizer
from lpplibs import MergeBackLPP
from lpplibs.AbbrGenerator import AbbrGenerator

t_start = time.clock()

pl_table = './lpplibs/CM_NormalLipids.xlsx'
fa_table = './lpplibs/FA_list.csv'
mod_table = './lpplibs/ModConfig.csv'

save_sdf = 'new_method_sdf_PA_max_1keto_1lessDB.sdf'
sdf_writer = Chem.SDWriter(save_sdf)
sdf_dct = {}

# pl_class_use_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']
pl_class_use_lst = ['PA']

parser = PLParser()
abbr_gen = AbbrGenerator()

pl_df = pd.read_excel(pl_table, sheetname=1)
fa_df = pd.read_csv(fa_table, index_col=0)

print(pl_df.head())
c_lst = []

fa_lpp_df_dct = {}

sum_theo_lpp_dct = {}
for (_idx, _row) in pl_df.iterrows():

    _pl_abbr = _row['phospholipids']

    _pl_elem_lst = parser.get_composition(_pl_abbr)
    print ('PL composition ==>', _pl_elem_lst)
    _pl_hg_abbr = _pl_elem_lst[0]

    # get smiles from abbr

    if _pl_hg_abbr in pl_class_use_lst:
        c_lst.append(_pl_abbr)

        # prepare output
        _pl_lpp_df = pd.DataFrame()

        print('Start oxidation of ==>', _pl_abbr)
        _pl_sn1_abbr = _pl_elem_lst[1]
        _pl_sn2_abbr = _pl_elem_lst[2]
        _pl_sn1_smiles = fa_df.loc[_pl_sn1_abbr, 'SMILES']
        _pl_sn2_smiles = fa_df.loc[_pl_sn2_abbr, 'SMILES']
        print('sn1 =>', _pl_sn1_smiles, '|| sn2 =>', _pl_sn2_smiles)

        # check if FA already oxidized to speed up
        if _pl_sn1_abbr in fa_lpp_df_dct.keys():
            sn1_mod_sum_df = fa_lpp_df_dct[_pl_sn1_abbr]
        else:
            sn1_link_dct = fa_link_filter(_pl_sn1_smiles)
            sn1_mod_sum_df = oxidizer(sn1_link_dct, mod_table)
            fa_lpp_df_dct[_pl_sn1_abbr] = sn1_mod_sum_df.copy()
            
        if _pl_sn2_abbr in fa_lpp_df_dct.keys():
            sn2_mod_sum_df = fa_lpp_df_dct[_pl_sn2_abbr]
        else:
            sn2_link_dct = fa_link_filter(_pl_sn2_smiles)
            sn2_mod_sum_df = oxidizer(sn2_link_dct, mod_table)
            fa_lpp_df_dct[_pl_sn2_abbr] = sn2_mod_sum_df.copy()

        for (_sn1_idx, _sn1_row) in sn1_mod_sum_df.iterrows():
            _sn1_mod_smiles = _sn1_row['FULL_SMILES']
            _sn1_abbr_str = _sn1_row['FA_ABBR']
            _sn1_typ_str = _sn1_row['FA_TYPE']
            # _sn1_json = _sn1_row['FA_JSON']
            # print(_sn1_row)
            # _sn1_abbr_str, _sn1_typ_str = abbr_gen.decode(_sn1_row['FA_CHECKER'])
            # print(_sn1_row)
            for (_sn2_idx, _sn2_row) in sn2_mod_sum_df.iterrows():
                _sn2_mod_smiles = _sn2_row['FULL_SMILES']
                _sn2_abbr_str = _sn2_row['FA_ABBR']
                _sn2_typ_str = _sn2_row['FA_TYPE']
                # _sn2_abbr_str, _sn2_typ_str = abbr_gen.decode(_sn2_row['FA_CHECKER'])
                # print(_sn2_row)

                # print('sn1 LPP =>', _sn1_abbr_str, ' | sn2 LPP =>', _sn2_abbr_str)

                _oap_ocp_lst = [_sn1_typ_str, _sn2_typ_str]
                _lpp_typ = ''.join(_oap_ocp_lst)

                # only export OAP & OCP
                if _lpp_typ not in ['LYSOLYSO', 'UNMODUNMOD']:
                    _lpp_smiles = MergeBackLPP.pl_lpp(_pl_hg_abbr, _sn2_mod_smiles, _sn1_mod_smiles)
                    _lpp_id_str = str(''.join([_pl_hg_abbr, '(', _sn1_abbr_str, '/', _sn2_abbr_str, ')']))

                    _lpp_sub_class_json = '{"SN1": "%s", "SN2": "%s"}' % (_sn1_typ_str, _sn2_typ_str)
                    # _lpp_name_str = _lpp_id_str.replace('/', '_')
                    # _lpp_name_str = _lpp_name_str.replace(':', '-')
                    # _lpp_name_str = _lpp_name_str.replace('@', 'at')
                    # _lpp_name_str = _lpp_name_str.replace('<', '[')
                    # _lpp_name_str = _lpp_name_str.replace('>', ']')
                    # print (_lpp_name_str)
                    _lpp_info_dct = {'LPP_ORIGIN': _pl_abbr, 'LPP_SMILES': _lpp_smiles, 'LPP_CLASS': _pl_hg_abbr,
                                     'SN1_SMILES': _sn1_mod_smiles, 'SN2_SMILES': _sn2_mod_smiles,
                                     'SN1_ABBR': _sn1_abbr_str, 'SN2_ABBR': _sn2_abbr_str,
                                     'SN1_JSON': _sn1_row['FA_JSON'], 'SN2_JSON': _sn2_row['FA_JSON'],
                                     'LM_ID': _lpp_id_str, 'SN_JSON': _lpp_sub_class_json}
                    # 'SN1_INFO': _sn1_row['FA_CHECKER'], 'SN2_INFO': _sn2_row['FA_CHECKER'],

                    _lpp_info_se = pd.Series(data=_lpp_info_dct)
                    _pl_lpp_df[_lpp_id_str] = _lpp_info_se

                    # check if same lpp generated already
                    # Currently use bulk settings
                    if _lpp_id_str in sdf_dct.keys():
                        _lpp_origin = sdf_dct[_lpp_id_str]['LPP_ORIGIN']
                        _lpp_origin_lst = _lpp_origin.split(',')
                        if _pl_abbr in _lpp_origin_lst:
                            pass
                        else:
                            _lpp_origin_lst.append(_pl_abbr)
                            sdf_dct[_lpp_id_str]['LPP_ORIGIN'] = ','.join(_lpp_origin_lst)
                    else:
                        sdf_dct[_lpp_id_str] = _lpp_info_dct.copy()

                    # clean memory by deleting these dicts and series
                    del _lpp_info_dct, _lpp_info_se

        # generate summary table
        _pl_lpp_df = _pl_lpp_df.transpose()
        print('==> %i of LPP generated !!' % _pl_lpp_df.shape[0])
        print('==> ==> Move to next lipid==> ')
        # print(_pl_lpp_df.head())
        sum_theo_lpp_dct[_pl_abbr] = _pl_lpp_df

        # create sdf
        # for (_lpp_i, _lpp_r) in _pl_lpp_df.iterrows():

sum_theo_lpp_pl = pd.Panel(data=sum_theo_lpp_dct)
print(sum_theo_lpp_pl.shape)

# write to sdf
print('==>Start to generate SDF ==>')
print('!! %i structures in total !!' % len(sdf_dct.keys()))
for _k_lpp in sdf_dct.keys():
    _lpp_dct = sdf_dct[_k_lpp]
    _lpp_smiles = str(_lpp_dct['LPP_SMILES'])
    # print(_lpp_smiles)
    _lpp_mol = Chem.MolFromSmiles(_lpp_smiles)
    AllChem.Compute2DCoords(_lpp_mol)
    _lpp_mol.SetProp('_Name', str(_lpp_dct['LM_ID']))
    _lpp_mass = Descriptors.MolWt(_lpp_mol)
    _lpp_exactmass = Descriptors.ExactMolWt(_lpp_mol)
    _lpp_formula = rdMolDescriptors.CalcMolFormula(_lpp_mol)
    _lpp_mol.SetProp('EXACT_MASS', '%.6f' % _lpp_exactmass)
    _lpp_mol.SetProp('NOMINAL_MASS', '%.3f' % _lpp_mass)
    _lpp_mol.SetProp('FORMULA', _lpp_formula)
    # _lpp_mol

    for _k in _lpp_dct.keys():
        _lpp_mol.SetProp(_k, str(_lpp_dct[_k]))

    sdf_writer.write(_lpp_mol)

sdf_writer.close()

t_spent = time.clock() - t_start
print('==>==>%i of LPP generated ==> ==> ' % len(sdf_dct.keys()))
print('==>==>==> %i of phospholipids processed in %.3fs ==> ==> ==> Finished !!!!!!' % (len(c_lst), t_spent))
