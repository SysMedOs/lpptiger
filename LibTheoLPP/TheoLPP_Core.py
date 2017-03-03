# -*- coding: utf-8 -*-
#
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPtiger.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#

from __future__ import print_function

import json
import time

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from LibTheoLPP import MSPcreator
from LibTheoLPP import LPPmerge
from LibTheoLPP import SDFsummary
from LibTheoLPP.AbbrGenerator import AbbrGenerator
from LibTheoLPP.TheoOxidation import fa_link_filter, oxidizer
from LibTheoLPP.ExactMassCalc import MZcalc
from LibTheoLPP.LPPparser import PLParser
from LibTheoLPP.TheoFrag import TheoFrag
from LibTheoLPP.FingerprintGenerator import FingerprintGen


def theolpp(usr_params):
    """
    param_dct = {'lipid_class': lipid_class, 'ox_level': ox_level,
                 'oap_mode': oap_mode, 'ocp_mode': ocp_mode,
                 'lyso_oap_mode': lyso_oap_mode, 'lyso_ocp_mode': lyso_ocp_mode,
                 'ox_max': ox_max, 'keto_max': keto_max, 'ooh_max': ooh_max, 'epoxy_max': epoxy_max,
                 'lipid_lst_path': lipid_lst_path, 'lipid_tab': lipid_tab,
                 'prostane_mode': prostane_mode, 'ox_prostane_mode': ox_prostane_mode,
                 'sdf_path': sdf_path, 'msp_mode': msp_mode, 'msp_path': msp_path,
                 'mod_lst_path': mod_lst_path, 'fa_lst_path': fa_lst_path, 'prostane_mod_path': prostane_mod_path,
                 'prostane_abbr_path': prostane_abbr_path, 'frag_pattern_path': frag_pattern_path}
    :param usr_params:
    :return:
    """

    t_start = time.clock()

    pl_table = usr_params['lipid_lst_path']
    fa_table = usr_params['fa_lst_path']
    mod_table = usr_params['mod_lst_path']
    isop_cfg = usr_params['prostane_mod_path']
    isopabbr_cfg = usr_params['prostane_abbr_path']
    # pl_class_use_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']
    pl_class = usr_params['lipid_class']
    pl_class_use_lst = [pl_class]
    ox_level = usr_params['ox_level']

    oap_mode = usr_params['oap_mode']
    ocp_mode = usr_params['ocp_mode']
    lyso_oap_mode = usr_params['lyso_oap_mode']
    lyso_ocp_mode = usr_params['lyso_ocp_mode']

    ox_max = usr_params['ox_max']
    keto_max = usr_params['keto_max']
    ooh_max = usr_params['ooh_max']
    epoxy_max = usr_params['epoxy_max']

    prostane_mode = usr_params['prostane_mode']
    prostane_ox_mode = usr_params['ox_prostane_mode']
    save_sdf = usr_params['sdf_path']
    save_spectra = usr_params['msp_mode']
    save_msp = usr_params['msp_path']
    score_xlsx = usr_params['frag_pattern_path']
    pl_fp_xlsx = usr_params['pl_hg_path']

    pl_df = pd.read_excel(pl_table, sheetname=usr_params['lipid_tab'])
    fa_df = pd.read_csv(fa_table, index_col=0)
    print(pl_df.head())

    # Select export species OAP, OCP, Lyso OAP, Lyso OCP
    ban_lst = ['LYSOLYSO']
    if oap_mode == 0:
        ban_lst.extend(['UNMODOAP', 'OAPUNMOD', 'OAPOAP'])
    if ocp_mode == 0:
        ban_lst.extend(['UNMODOCP', 'OCPUNMOD', 'OCPOCP'])
    if lyso_oap_mode == 0:
        ban_lst.extend(['LYSOOAP', 'OAPLYSO'])
    if lyso_ocp_mode == 0:
        ban_lst.extend(['LYSOOCP', 'OCPLYSO'])
    if ox_level == 1:
        ban_lst.extend(['OAPOAP', 'OCPOCP', 'OAPOCP', 'OCPOAP', 'OAPUNMOD', 'OCPUNMOD'])

    ox_param_dct = {'MAX_MOD': ox_max, 'MAX_KETO': keto_max, 'MAX_OOH': ooh_max, 'MAX_EPOXY': epoxy_max}

    sdf_writer = Chem.SDWriter(open(save_sdf, mode='w'))
    if save_spectra == 1 and len(save_msp) > 0:
        msp_obj = open(save_msp, mode='w')
    else:
        msp_obj = None
    sdf_dct = {}

    parser = PLParser()
    abbr_gen = AbbrGenerator()

    frag_gen = TheoFrag(pl_class, score_xlsx)
    fingerprint_gen = FingerprintGen(pl_fp_xlsx)

    c_lst = []

    fa_lpp_df_dct = {}

    sum_theo_lpp_dct = {}
    for (_idx, _row) in pl_df.iterrows():

        _pl_abbr = str(_row['phospholipids'])

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
                sn1_mod_sum_df = oxidizer(sn1_link_dct, mod_table, isop_cfg, isopabbr_cfg,
                                          ox_level, ox_param_dct, prostane_mode, prostane_ox_mode)
                fa_lpp_df_dct[_pl_sn1_abbr] = sn1_mod_sum_df.copy()

            if _pl_sn2_abbr in fa_lpp_df_dct.keys():
                sn2_mod_sum_df = fa_lpp_df_dct[_pl_sn2_abbr]
            else:
                sn2_link_dct = fa_link_filter(_pl_sn2_smiles)
                sn2_mod_sum_df = oxidizer(sn2_link_dct, mod_table, isop_cfg, isopabbr_cfg,
                                          ox_level, ox_param_dct, prostane_mode, prostane_ox_mode)
                fa_lpp_df_dct[_pl_sn2_abbr] = sn2_mod_sum_df.copy()

            for (_sn1_idx, _sn1_row) in sn1_mod_sum_df.iterrows():
                _sn1_mod_smiles = _sn1_row['FULL_SMILES']
                _sn1_abbr_str = _sn1_row['FA_ABBR']
                _sn1_typ_str = _sn1_row['FA_TYPE']
                _sn1_formula_str = _sn1_row['FA_FORMULA']

                for (_sn2_idx, _sn2_row) in sn2_mod_sum_df.iterrows():
                    _sn2_mod_smiles = _sn2_row['FULL_SMILES']
                    _sn2_abbr_str = _sn2_row['FA_ABBR']
                    _sn2_typ_str = _sn2_row['FA_TYPE']
                    _sn2_formula_str = _sn2_row['FA_FORMULA']

                    _oap_ocp_lst = [_sn1_typ_str, _sn2_typ_str]
                    _lpp_typ = ''.join(_oap_ocp_lst)

                    if _lpp_typ not in ban_lst:
                        _lpp_smiles = LPPmerge.pl_lpp(_pl_hg_abbr, sn1=_sn1_mod_smiles, sn2=_sn2_mod_smiles)
                        _lpp_id_str = str(''.join([_pl_hg_abbr, '(', _sn1_abbr_str, '/', _sn2_abbr_str, ')']))

                        _lpp_sub_class_json = '{"SN1": "%s", "SN2": "%s"}' % (_sn1_typ_str, _sn2_typ_str)

                        _lpp_info_dct = {'LPP_ORIGIN': _pl_abbr, 'LPP_SMILES': _lpp_smiles, 'LPP_CLASS': _pl_hg_abbr,
                                         'SN1_SMILES': _sn1_mod_smiles, 'SN2_SMILES': _sn2_mod_smiles,
                                         'SN1_ABBR': _sn1_abbr_str, 'SN2_ABBR': _sn2_abbr_str,
                                         'SN1_JSON': _sn1_row['FA_JSON'], 'SN2_JSON': _sn2_row['FA_JSON'],
                                         'SN1_FORMULA': _sn1_formula_str, 'SN2_FORMULA': _sn2_formula_str,
                                         'LM_ID': _lpp_id_str, 'SN_JSON': _lpp_sub_class_json}
                        if save_spectra == 1:
                            _lpp_info_dct['MSP_JSON'] = frag_gen.calc_frags(_lpp_info_dct)
                            print('_lpp_info_dct[MSP_JSON]')
                            print(_lpp_info_dct['MSP_JSON'])

                        # print(_lpp_info_dct)
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
    print('==>Start to generate SDF ==> MSP mode = %i' % save_spectra)
    print('!! %i structures in total !!' % len(sdf_dct.keys()))

    mzcalc = MZcalc()

    if save_spectra == 1:
        for _k_lpp in sdf_dct.keys():
            _lpp_dct = sdf_dct[_k_lpp]
            if len(json.loads(_lpp_dct['MSP_JSON']).keys()) > 0:
                _lpp_smiles = str(_lpp_dct['LPP_SMILES'])
                # print(_lpp_smiles)
                _lpp_mol = Chem.MolFromSmiles(_lpp_smiles)
                AllChem.Compute2DCoords(_lpp_mol)
                _lpp_mol.SetProp('_Name', str(_lpp_dct['LM_ID']))
                _lpp_mass = Descriptors.MolWt(_lpp_mol)
                _lpp_exactmass = rdMolDescriptors.CalcExactMolWt(_lpp_mol)
                _lpp_formula = rdMolDescriptors.CalcMolFormula(_lpp_mol)
                _lpp_mol.SetProp('EXACT_MASS', '%.6f' % _lpp_exactmass)
                _lpp_mol.SetProp('NOMINAL_MASS', '%.3f' % _lpp_mass)
                _lpp_mol.SetProp('FORMULA', _lpp_formula)
                _lpp_sn2_smi = _lpp_dct['SN2_SMILES']

                if str(_lpp_dct['LPP_CLASS']) == 'PC' and _lpp_sn2_smi[-9:] != r'C(O)=O)=O':
                    _lpp_neg_precursor_elem = mzcalc.get_elements(_lpp_formula)
                    _lpp_neg_precursor_formula = mzcalc.get_formula(_lpp_neg_precursor_elem, charge='[M+HCOO]-')
                    _lpp_neg_precursor_mz = mzcalc.get_mono_mz(_lpp_formula, charge='[M+HCOO]-')
                    _lpp_neg_precursor_info = '{"[M+HCOO]-": ["%s", %f]}' % (_lpp_neg_precursor_formula[0],
                                                                             _lpp_neg_precursor_mz)

                else:
                    _lpp_neg_precursor_elem = mzcalc.get_elements(_lpp_formula)
                    _lpp_neg_precursor_formula = mzcalc.get_formula(_lpp_neg_precursor_elem, charge='[M-H]-')
                    _lpp_neg_precursor_mz = mzcalc.get_mono_mz(_lpp_formula, charge='[M-H]-')
                    _lpp_neg_precursor_info = '{"[M-H]-": ["%s", %f]}' % (_lpp_neg_precursor_formula[0],
                                                                          _lpp_neg_precursor_mz)

                _lpp_dct['PRECURSOR_JSON'] = _lpp_neg_precursor_info
                _lpp_mol.SetProp('PRECURSOR_JSON', _lpp_neg_precursor_info)
                _lpp_dct['EXACT_MASS'] = _lpp_exactmass
                fp_mz_lst = fingerprint_gen.get_fingerprint(_lpp_dct)
                _lpp_dct['FINGERPRINT'] = fp_mz_lst
                _lpp_mol.SetProp('FINGERPRINT', json.dumps(fp_mz_lst))

                for _k in _lpp_dct.keys():
                    _lpp_mol.SetProp(_k, str(_lpp_dct[_k]))

                sdf_writer.write(_lpp_mol)
                if save_spectra == 1 and len(save_msp) > 0:
                    MSPcreator.to_msp(msp_obj, _lpp_dct)

    elif save_spectra == 0:
        for _k_lpp in sdf_dct.keys():
            _lpp_dct = sdf_dct[_k_lpp]
            _lpp_smiles = str(_lpp_dct['LPP_SMILES'])
            print(_lpp_smiles)
            _lpp_mol = Chem.MolFromSmiles(_lpp_smiles)
            AllChem.Compute2DCoords(_lpp_mol)
            _lpp_mol.SetProp('_Name', str(_lpp_dct['LM_ID']))
            _lpp_mass = Descriptors.MolWt(_lpp_mol)
            _lpp_exactmass = rdMolDescriptors.CalcExactMolWt(_lpp_mol)
            _lpp_formula = rdMolDescriptors.CalcMolFormula(_lpp_mol)
            _lpp_mol.SetProp('EXACT_MASS', '%.6f' % _lpp_exactmass)
            _lpp_mol.SetProp('NOMINAL_MASS', '%.3f' % _lpp_mass)
            _lpp_mol.SetProp('FORMULA', _lpp_formula)
            _lpp_sn2_smi = _lpp_dct['SN2_SMILES']

            if str(_lpp_dct['LPP_CLASS']) == 'PC' and _lpp_sn2_smi[-9:] != r'C(O)=O)=O':
                _lpp_neg_precursor_elem = mzcalc.get_elements(_lpp_formula)
                _lpp_neg_precursor_formula = mzcalc.get_formula(_lpp_neg_precursor_elem, charge='[M+HCOO]-')
                _lpp_neg_precursor_mz = mzcalc.get_mono_mz(_lpp_formula, charge='[M+HCOO]-')
                _lpp_neg_precursor_info = '{"[M+HCOO]-": ["%s", %f]}' % (_lpp_neg_precursor_formula[0],
                                                                         _lpp_neg_precursor_mz)

            else:
                _lpp_neg_precursor_elem = mzcalc.get_elements(_lpp_formula)
                _lpp_neg_precursor_formula = mzcalc.get_formula(_lpp_neg_precursor_elem, charge='[M-H]-')
                _lpp_neg_precursor_mz = mzcalc.get_mono_mz(_lpp_formula, charge='[M-H]-')
                _lpp_neg_precursor_info = '{"[M-H]-": ["%s", %f]}' % (_lpp_neg_precursor_formula[0],
                                                                      _lpp_neg_precursor_mz)

            _lpp_dct['PRECURSOR_JSON'] = _lpp_neg_precursor_info
            _lpp_mol.SetProp('PRECURSOR_JSON', _lpp_neg_precursor_info)

            for _k in _lpp_dct.keys():
                _lpp_mol.SetProp(_k, str(_lpp_dct[_k]))

            sdf_writer.write(_lpp_mol)

    sdf_writer.close()
    if save_spectra == 1 and len(save_msp) > 0:
        msp_obj.close()

    SDFsummary.sdf2xlsx(save_sdf, str(save_sdf)[:-4] + '.xlsx')
    if save_msp == 1:
        SDFsummary.sdf2sum_fa(save_sdf, str(save_sdf)[:-4] + '_FA_SUM.xlsx')

    t_spent = time.clock() - t_start
    info_updater_1 = '=>%i of LPP generated ==> ' % len(sdf_dct.keys())
    info_updater_2 = '=>==> %i of phospholipids processed in %.3fs ==> ==> Finished !!!!!!' % (len(c_lst), t_spent)

    return info_updater_1, info_updater_2

