# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import time
import json

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


class SDFparserPCA(object):

    def __init__(self, sdf_path_lst):
        self.sdf_path_lst = sdf_path_lst

    def export_for_pca(self, save_path):

        sdf_info_dct = {}

        for _sdf in self.sdf_path_lst:
            _sdf_obj = Chem.SDMolSupplier(_sdf)

            for _mol in _sdf_obj:

                _id = _mol.GetProp('LM_ID')
                _lpp_class = _mol.GetProp('LPP_CLASS')
                _sn_json = _mol.GetProp('SN_JSON')
                _sn1_json = _mol.GetProp('SN1_JSON')
                _sn2_json = _mol.GetProp('SN2_JSON')

                _sn_dct = json.loads(_sn_json)
                _sn1_dct = json.loads(_sn1_json)
                _sn2_dct = json.loads(_sn2_json)

                _mol_dct = {'ID': _id, 'LPP_CLASS': _lpp_class}

                for _sn1_key in _sn1_dct.keys():
                    _mol_dct['SN1_%s' % _sn1_key] = _sn1_dct[_sn1_key]
                for _sn2_key in _sn2_dct.keys():
                    _mol_dct['SN2_%s' % _sn2_key] = _sn2_dct[_sn2_key]

                for _sn_key in _sn_dct.keys():

                    if _sn_dct[_sn_key] == 'OAP':
                        _mol_dct['%s_OAP' % _sn_key] = 1
                        _mol_dct['%s_OCP' % _sn_key] = 0
                        _mol_dct['%s_UNMOD' % _sn_key] = 0
                        _mol_dct['%s_LYSO' % _sn_key] = 0
                    if _sn_dct[_sn_key] == 'OCP':
                        _mol_dct['%s_OAP' % _sn_key] = 0
                        _mol_dct['%s_OCP' % _sn_key] = 1
                        _mol_dct['%s_UNMOD' % _sn_key] = 0
                        _mol_dct['%s_LYSO' % _sn_key] = 0
                    if _sn_dct[_sn_key] == 'UNMOD':
                        _mol_dct['%s_OAP' % _sn_key] = 0
                        _mol_dct['%s_OCP' % _sn_key] = 0
                        _mol_dct['%s_UNMOD' % _sn_key] = 1
                        _mol_dct['%s_LYSO' % _sn_key] = 0
                    if _sn_dct[_sn_key] == 'LYSO':
                        _mol_dct['%s_OAP' % _sn_key] = 0
                        _mol_dct['%s_OCP' % _sn_key] = 0
                        _mol_dct['%s_UNMOD' % _sn_key] = 0
                        _mol_dct['%s_LYSO' % _sn_key] = 1

                sdf_info_dct[_id] = _mol_dct
                del _mol_dct

        sdf_info_df = pd.DataFrame(data=sdf_info_dct)
        sdf_info_df = sdf_info_df.transpose()
        sdf_info_df.to_excel(save_path)

        print('Finished!')

if __name__ == '__main__':

    pa_sdf_file = r'D:\theolpp\PA_CM_max_1keto_1lessDB_FRAG.sdf'
    pc_sdf_file = r'D:\theolpp\PC_CM_max_1keto_1lessDB_FRAG.sdf'
    pe_sdf_file = r'D:\theolpp\PE_CM_max_1keto_1lessDB_FRAG.sdf'
    pg_sdf_file = r'D:\theolpp\PG_CM_max_1keto_1lessDB_FRAG.sdf'
    ps_sdf_file = r'D:\theolpp\PS_CM_max_1keto_1lessDB_FRAG.sdf'
    save_as = r'D:\theolpp\PL_CM_max_1keto_1lessDB_FRAG_pca.xlsx'

    pl_sdf_lst = [pa_sdf_file, pc_sdf_file, pe_sdf_file, pg_sdf_file, ps_sdf_file]

    sdf_extractor = SDFparserPCA(sdf_path_lst=pl_sdf_lst)
    sdf_extractor.export_for_pca(save_path=save_as)
