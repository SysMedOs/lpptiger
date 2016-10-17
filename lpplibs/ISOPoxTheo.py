# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re
import json

import numpy as np
import pandas as pd

from AbbrGenerator import AbbrGenerator


class IsoProstanOx(object):

    def __init__(self, db_info_dct, isop_mod_cfg):
        self.db_info_dct = db_info_dct
        self.isop_mod_df = pd.read_csv(isop_mod_cfg, header=0, index_col=0)

    @staticmethod
    def get_count(isop_smi):

        c_count = isop_smi.count('C')


    @staticmethod
    def get_isop_abbr(isop_info_dct, usr_c, usr_db):

        # usr_code = 'P-18:0[0xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'

        abbr_formater = ('{RING_TYPE}{DB_COUNT}-{C_COUNT}:{DB_COUNT}[{DB_COUNT}xDB,{OH_COUNT}xOH,{KETO_COUNT}xKETO]'
                         '<CHO@C{CHO_COUNT},COOH@{COOH_COUNT}>{left}OAP:{OAP_COUNT},OCP:{OCP_COUNT}{right}'
                         .format(RING_TYPE=isop_info_dct['RING_TYPE'],
                                 C_COUNT=usr_c,
                                 DB_COUNT=usr_db,
                                 OH_COUNT=isop_info_dct['OH'],
                                 KETO_COUNT=isop_info_dct['KETO'],
                                 CHO_COUNT=isop_info_dct['CHO'],
                                 COOH_COUNT=isop_info_dct['COOH'],
                                 left='{',
                                 OAP_COUNT=isop_info_dct['OAP'],
                                 OCP_COUNT=isop_info_dct['OCP'],
                                 right='}'
                                 )
                         )

        print(abbr_formater)
        return abbr_formater

    def get_isop_lpp(self):

        # use only less C=C
        usr_db_num = self.db_info_dct['DB_count']
        usr_isop_mod_df = self.isop_mod_df[self.isop_mod_df['DB_NUM_REQUIRED'] <= usr_db_num]

        usr_db_main_smi = self.db_info_dct['DB_main_part']
        usr_db_pre_smi = self.db_info_dct['DB_pre_part']
        usr_db_post_smi = self.db_info_dct['DB_post_part']

        usr_lpp_df = pd.DataFrame()
        usr_lpp_dct = {}
        usr_lpp_lst = []

        # define the empty possibilities in configure .csv
        no_rev_lst = ['', 'NA', 'NULL', r'N/A', 'NAN', 'na', 'null', r'n/a', 'nan']

        for _abbr, _ring_info_se in usr_isop_mod_df.iterrows():

            _ring_info_dct = _ring_info_se.to_dict()

            print(_ring_info_dct['ORIGIN_SMILES'])
            print(usr_db_main_smi)
            _pattern_rgx = re.compile(_ring_info_dct['ORIGIN_SMILES'])

            _db_required_num = _ring_info_dct['DB_NUM_REQUIRED']
            _remain_db_num = usr_db_num - _db_required_num
            _pre_ring_db_count_lst = range(0, _remain_db_num + 1)
            _post_ring_db_count_lst = sorted(_pre_ring_db_count_lst, reverse=True)
            _pre_post_db_count_lst = zip(_pre_ring_db_count_lst, _post_ring_db_count_lst)
            print(_pre_post_db_count_lst)
            for _db_seg in _pre_post_db_count_lst:

                _rebuild_smi = '/C=C\\C' * _db_seg[0] + _ring_info_dct['ORIGIN_SMILES'] + '/C=C\\C' * _db_seg[1]

                if _rebuild_smi == usr_db_main_smi:

                    # print('rebuild check pass!')
                    if _ring_info_dct['OCP'] == 0 and _ring_info_dct['OAP'] == 1:
                        _ring_lpp_smi = (usr_db_pre_smi + '/C=C\\C' * _db_seg[0] +
                                         _ring_info_dct['MAIN_SMILES'] + '/C=C\\C' * _db_seg[1] + usr_db_post_smi)
                        _c_count = _ring_lpp_smi.count('C')
                        print(_ring_info_dct)
                        _abbr = self.get_isop_abbr(_ring_info_dct, usr_c=_c_count, usr_db=_remain_db_num)
                        usr_lpp_lst.append(_ring_lpp_smi)
                    elif _ring_info_dct['OCP'] == 1 and _ring_info_dct['OAP'] == 0:
                        _ring_lpp_smi = (usr_db_pre_smi + '/C=C\\C' * _db_seg[0] +
                                         _ring_info_dct['MAIN_SMILES'] + usr_db_post_smi)
                        usr_lpp_lst.append(_ring_lpp_smi)

                    _rev_smi = _ring_info_dct['REVERSE_SMILES']
                    if isinstance(_rev_smi, str):
                        # blank 'REVERSE_SMILES' gives np.nan which is float
                        if _rev_smi not in no_rev_lst:
                            if _ring_info_dct['OCP'] == 0 and _ring_info_dct['OAP'] == 1:
                                _rev_ring_lpp_smi = (usr_db_pre_smi + '/C=C\\C' * _db_seg[0] +
                                                     _ring_info_dct['REVERSE_SMILES']
                                                     + '/C=C\\C' * _db_seg[1]+ usr_db_post_smi)

                                usr_lpp_lst.append(_rev_ring_lpp_smi)
                            elif _ring_info_dct['OCP'] == 1 and _ring_info_dct['OAP'] == 0:
                                _rev_ring_lpp_smi = (usr_db_pre_smi + '/C=C\\C' * _db_seg[0] +
                                                     _ring_info_dct['REVERSE_SMILES'] + usr_db_post_smi)
                                usr_lpp_lst.append(_rev_ring_lpp_smi)
                        else:
                            pass
                    else:
                        pass

        # ox_isop_lpp_lst = [usr_db_pre_smi + _lpp_smi + usr_db_post_smi for _lpp_smi in usr_lpp_lst]

        print(usr_lpp_lst)

        return usr_lpp_lst


if __name__ == '__main__':

    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem

    isop_cfg = r'IsoP_ModConfig.csv'

    ara_smi = r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O'
    epa_smi = r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CC)=O'

    linolenic_db_info_dct = {'DB_count': 3,
                             'DB_pre_part': r'OC(CCCCCCC',
                             'DB_main_part': r'/C=C\C/C=C\C/C=C\C',
                             'DB_post_part': r'C)=O',
                             'DB_full_fa': r'OC(CCCCCCC/C=C\C/C=C\C/C=C\CC)=O',
                             'DB_LINK_type': '', 'DB_C_count': 18, 'DB_start': 9, 'DB_omega': 3}

    ara_db_info_dct = {'DB_count': 4,
                       'DB_pre_part': r'OC(CCC',
                       'DB_main_part': r'/C=C\C/C=C\C/C=C\C/C=C\C',
                       'DB_post_part': r'CCCCC)=O',
                       'DB_full_fa': r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O',
                       'DB_LINK_type': '', 'DB_C_count': 20, 'DB_start': 5, 'DB_omega': 6}

    epa_db_info_dct = {'DB_count': 5,
                       'DB_pre_part': r'OC(CCC',
                       'DB_main_part': r'/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C',
                       'DB_post_part': r'C)=O',
                       'DB_full_fa': r'OC(CCC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CC)=O',
                       'DB_LINK_type': '', 'DB_C_count': 20, 'DB_start': 5, 'DB_omega': 3}

    save_img = 'test_epa_isop_lpp.png'
    save_sdf = 'test_epa_IsoP.sdf'

    ox_isop = IsoProstanOx(epa_db_info_dct, isop_cfg)

    got_lpp_lst = ox_isop.get_isop_lpp()

    lpp_mol_lst = []

    sdf_writer = Chem.SDWriter(save_sdf)

    for _isop_lpp in got_lpp_lst:
        try:
            _mol = Chem.MolFromSmiles(_isop_lpp)

            AllChem.Compute2DCoords(_mol)
            sdf_writer.write(_mol)
            lpp_mol_lst.append(_mol)
        except:
            error_rgx = re.compile(r'//')
            error_checker = error_rgx.sub(r'', _isop_lpp)
            # print(error_checker)
            try:
                _mol = Chem.MolFromSmiles(error_checker)
                AllChem.Compute2DCoords(_mol)
                lpp_mol_lst.append(_mol)
                print('corrected and pass')
            except:
                print('not working for --> ', _isop_lpp)

            # pass

    print(len(lpp_mol_lst))
    img = Draw.MolsToGridImage(lpp_mol_lst, molsPerRow=6, subImgSize=(200, 200))
    # img.show()
    img.save(save_img)
    print('saved!')
