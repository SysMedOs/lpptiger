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
from lpplibs.SNMainFrag import SNMainFrag
from lpplibs.AbbrGenerator import AbbrGenerator
from lpplibs import MSPcreator
from lpplibs import SDFsummary


def theolpp(usr_params):
    """
    param_dct = {'lipid_class': lipid_class, 'ox_level': ox_level, 'ox_max': ox_max,
                 'lipid_lst_path': lipid_lst_path, 'lipid_tab': lipid_tab,
                 'prostane_mode': prostane_mode, 'ox_prostane_mode': ox_prostane_mode,
                 'sdf_path': sdf_path, 'msp_mode': msp_mode, 'msp_path': msp_path,
                 'mod_lst_path': mod_lst_path, 'fa_lst_path': fa_lst_path, 'prostane_mod_path': prostane_mod_path,
                 'prostane_abbr_path': prostane_abbr_path, 'frag_pattern_path': frag_pattern_path}
    :param params:
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
    ox_max = usr_params['ox_max']
    save_sdf = usr_params['sdf_path']
    save_spectra = usr_params['msp_mode']
    save_msp = usr_params['msp_path']
    score_xlsx = usr_params['frag_pattern_path']

    pl_df = pd.read_excel(pl_table, sheetname=usr_params['lipid_tab'])
    fa_df = pd.read_csv(fa_table, index_col=0)

    # pl_table = r'./lpplibs/DW_PL.xlsx'
    # fa_table = r'./lpplibs/FA_list.csv'
    # mod_table = r'./lpplibs/ModConfig.csv'
    #
    # isop_cfg = r'./lpplibs/IsoP_ModConfig.csv'
    # isopabbr_cfg = r'./lpplibs/IsoP_AbbrConfig.csv'
    #
    # # pl_class_use_lst = ['PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PS']
    # pl_class_use_lst = ['PC']
    # pl_class = pl_class_use_lst[0]
    # ox_level = 1
    # save_spectra = 1
    #
    # save_sdf = '%s_DW_lv1.sdf' % ''.join(pl_class_use_lst)
    # save_msp = '%s_DW_lv1.msp' % ''.join(pl_class_use_lst)
    #
    # score_xlsx = './TheoFragPatterns_csv/ion_scores_df.xlsx'
    #
    # pl_df = pd.read_excel(pl_table, sheetname=0)
    # fa_df = pd.read_csv(fa_table, index_col=0)
    # max_o = 3

    print(pl_df.head())

    sdf_writer = Chem.SDWriter(open(save_sdf, mode='w'))
    msp_obj = open(save_msp, mode='w')
    sdf_dct = {}

    parser = PLParser()
    abbr_gen = AbbrGenerator()

    frag_gen = SNMainFrag(pl_class, score_xlsx)


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
                sn1_mod_sum_df = oxidizer(sn1_link_dct, mod_table, isop_cfg, isopabbr_cfg, ox_level)
                fa_lpp_df_dct[_pl_sn1_abbr] = sn1_mod_sum_df.copy()

            if _pl_sn2_abbr in fa_lpp_df_dct.keys():
                sn2_mod_sum_df = fa_lpp_df_dct[_pl_sn2_abbr]
            else:
                sn2_link_dct = fa_link_filter(_pl_sn2_smiles)
                sn2_mod_sum_df = oxidizer(sn2_link_dct, mod_table, isop_cfg, isopabbr_cfg, ox_level)
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
                    # if _lpp_typ not in ['LYSOLYSO', 'UNMODUNMOD']:
                    if _lpp_typ not in ['LYSOLYSO']:
                        _lpp_smiles = MergeBackLPP.pl_lpp(_pl_hg_abbr, sn1=_sn1_mod_smiles, sn2=_sn2_mod_smiles)
                        _lpp_id_str = str(''.join([_pl_hg_abbr, '(', _sn1_abbr_str, '/', _sn2_abbr_str, ')']))

                        _lpp_sub_class_json = '{"SN1": "%s", "SN2": "%s"}' % (_sn1_typ_str, _sn2_typ_str)

                        # _lpp_sn1_frag_lst = json.loads(_sn1_row['FRAG_SMILES'])
                        # if _lpp_sn1_frag_lst != ['']:
                        #     _lpp_sn1_frag_lst.append(_sn1_row['FULL_SMILES'])
                        # else:
                        #     _lpp_sn1_frag_lst = [_sn1_row['FULL_SMILES']]
                        #
                        # _lpp_sn2_frag_lst = json.loads(_sn2_row['FRAG_SMILES'])
                        # if _lpp_sn2_frag_lst != ['']:
                        #     _lpp_sn2_frag_lst.append(_sn2_row['FULL_SMILES'])
                        # else:
                        #     _lpp_sn2_frag_lst = [_sn2_row['FULL_SMILES']]
                        #
                        # _lpp_frag_lst = []
                        #
                        # for _sn1_frag in _lpp_sn1_frag_lst:
                        #     for _sn2_frag in _lpp_sn2_frag_lst:
                        #         _lpp_frag_lst.append(MergeBackLPP.pl_lpp(_pl_hg_abbr,
                        #                                                  sn1=_sn1_frag, sn2=_sn2_frag))
                        # _lpp_frag_json = json.dumps(_lpp_frag_lst)

                        # 'LPP_FRAG': _lpp_frag_json,
                        # 'SN1_FRAGS': _sn1_row['FRAG_SMILES'], 'SN2_FRAGS': _sn2_row['FRAG_SMILES'],
                        _lpp_info_dct = {'LPP_ORIGIN': _pl_abbr, 'LPP_SMILES': _lpp_smiles, 'LPP_CLASS': _pl_hg_abbr,
                                         'SN1_SMILES': _sn1_mod_smiles, 'SN2_SMILES': _sn2_mod_smiles,
                                         'SN1_ABBR': _sn1_abbr_str, 'SN2_ABBR': _sn2_abbr_str,
                                         'SN1_JSON': _sn1_row['FA_JSON'], 'SN2_JSON': _sn2_row['FA_JSON'],
                                         'LM_ID': _lpp_id_str, 'SN_JSON': _lpp_sub_class_json}
                        if save_spectra == 1:
                            _lpp_info_dct['MSP_JSON'] = frag_gen.calc_frags(_lpp_info_dct)

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
    print('==>Start to generate SDF ==>')
    print('!! %i structures in total !!' % len(sdf_dct.keys()))
    for _k_lpp in sdf_dct.keys():

        _lpp_dct = sdf_dct[_k_lpp]

        if save_spectra == 1 and len(json.loads(_lpp_dct['MSP_JSON']).keys()) > 0:
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

            if str(_lpp_dct['LPP_CLASS']) == 'PC':
                _lpp_neg_precursor_mz = frag_gen.formula_to_mz(_lpp_formula, charge='[M+FA-H]-')
                _lpp_neg_precursor_info = '{"[M+FA-H]-": ["%s", %f]}' % (_lpp_formula, _lpp_neg_precursor_mz[0])

            else:
                _lpp_neg_precursor_mz = frag_gen.formula_to_mz(_lpp_formula, charge='[M-H]-')
                _lpp_neg_precursor_info = '{"[M-H]-": ["%s", %f]}' % (_lpp_formula, _lpp_neg_precursor_mz[0])

            _lpp_dct['PRECURSOR_JSON'] = _lpp_neg_precursor_info
            _lpp_mol.SetProp('PRECURSOR_JSON', _lpp_neg_precursor_info)

            for _k in _lpp_dct.keys():
                _lpp_mol.SetProp(_k, str(_lpp_dct[_k]))

            sdf_writer.write(_lpp_mol)
            MSPcreator.to_msp(msp_obj, _lpp_dct)

        elif save_spectra == 0:
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

            if str(_lpp_dct['LPP_CLASS']) == 'PC':
                _lpp_neg_precursor_mz = frag_gen.formula_to_mz(_lpp_formula, charge='[M+FA-H]-')
                _lpp_neg_precursor_info = '{"[M+FA-H]-": ["%s", %f]}' % (_lpp_formula, _lpp_neg_precursor_mz[0])

            else:
                _lpp_neg_precursor_mz = frag_gen.formula_to_mz(_lpp_formula, charge='[M-H]-')
                _lpp_neg_precursor_info = '{"[M-H]-": ["%s", %f]}' % (_lpp_formula, _lpp_neg_precursor_mz[0])

            _lpp_dct['PRECURSOR_JSON'] = _lpp_neg_precursor_info
            _lpp_mol.SetProp('PRECURSOR_JSON', _lpp_neg_precursor_info)

            for _k in _lpp_dct.keys():
                _lpp_mol.SetProp(_k, str(_lpp_dct[_k]))

            sdf_writer.write(_lpp_mol)

    sdf_writer.close()
    msp_obj.close()

    SDFsummary.sdf2xlsx(save_sdf, str(save_sdf)[:-4] + '.xlsx')

    t_spent = time.clock() - t_start
    print('==>==>%i of LPP generated ==> ==> ' % len(sdf_dct.keys()))
    print('==>==>==> %i of phospholipids processed in %.3fs ==> ==> ==> Finished !!!!!!' % (len(c_lst), t_spent))

