# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPsmi.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re
import json

import numpy as np
import pandas as pd

from AbbrGenerator import AbbrGenerator


class IsoProstanOx(object):
    """
    Create prostanes from PUFA
    """
    def __init__(self, db_info_dct, isop_mod_cfg, isop_abbr_cfg):
        """

        :param dict db_info_dct:
        :param str isop_mod_cfg:
        :param str isop_abbr_cfg:
        """
        self.db_info_dct = db_info_dct
        self.isop_mod_df = pd.read_csv(isop_mod_cfg, header=0, index_col=0)
        self.isop_abbr_df = pd.read_csv(isop_abbr_cfg, header=0)
        print(self.isop_abbr_df)

        # for _abbr_idx, _abbr_r in self.isop_abbr_df.iterrows():
        num_c_lst = [int(nc) for nc in self.isop_abbr_df['NUM_C'].tolist()]
        num_db_lst = [int(ndb) for ndb in self.isop_abbr_df['NUM_DB'].tolist()]
        class_abbr_lst = self.isop_abbr_df['CLASS_ABBR'].tolist()
        fa_pr_lst = zip(num_c_lst, num_db_lst)

        self.isop_abbr_dct = {}
        for _fa in fa_pr_lst:
            _fa_idx = fa_pr_lst.index(_fa)
            self.isop_abbr_dct[_fa] = class_abbr_lst[_fa_idx]

        # print('isop_abbr_dct', self.isop_abbr_dct)

    @staticmethod
    def get_isop_series(isop_info_dct, usr_pre_smi, usr_pre_db_count, reverse=False):

        pre_c_count = usr_pre_smi.count('C')

        if reverse is False:
            isop_series = pre_c_count + 1 + 3 * usr_pre_db_count + isop_info_dct['MAIN_SERIES']
            isop_series = str(isop_series)
        elif reverse is True:
            isop_series = pre_c_count + 1 + 3 * usr_pre_db_count + isop_info_dct['REVERSE_SERIES']
            isop_series = str(isop_series)
        else:
            isop_series = ''

        return isop_series

    def get_isop_checker(self, isop_info_dct, usr_c, usr_db, usr_series):

        # usr_code = 'P-18:0[0xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'
        if (usr_c, self.db_info_dct['DB_count']) in self.isop_abbr_dct.keys():
            isop_abbr = self.isop_abbr_dct[(usr_c, self.db_info_dct['DB_count'])]

            ring_type = isop_info_dct['RING_TYPE']
            class_type = isop_info_dct['CLASS_TYPE']

            # for [ABCDEFGHIJ]-IsoP and [DE]-IsoK
            if class_type in ['P', 'K']:
                abbr_formatter = ('{ISOP_SERIES}-{RING_TYPE}{SERIES_DB_COUNT}-{ISOP_ABBR}{ISOX}-'
                                  '{C_COUNT}:{DB_COUNT}[{DB_COUNT}xDB,'
                                  '{OH_COUNT}xOH,{KETO_COUNT}xKETO,{OOH_COUNT}xOOH]'
                                  '<CHO@C{CHO_COUNT},COOH@{COOH_COUNT}>{left}OAP:{OAP_COUNT},OCP:{OCP_COUNT}{right}'
                                  .format(ISOP_SERIES=usr_series,
                                          RING_TYPE=ring_type,
                                          SERIES_DB_COUNT=usr_db+isop_info_dct['DB'],
                                          ISOP_ABBR=isop_abbr,
                                          ISOX=class_type,
                                          C_COUNT=usr_c,
                                          DB_COUNT=usr_db+isop_info_dct['DB'],
                                          OH_COUNT=isop_info_dct['OH'],
                                          KETO_COUNT=isop_info_dct['KETO'],
                                          OOH_COUNT=isop_info_dct['OOH'],
                                          CHO_COUNT=isop_info_dct['CHO'],
                                          COOH_COUNT=isop_info_dct['COOH'],
                                          left='{',
                                          OAP_COUNT=isop_info_dct['OAP'],
                                          OCP_COUNT=isop_info_dct['OCP'],
                                          right='}'
                                          )
                                  )
                abbr_short = ('{RING_TYPE}{SERIES_DB_COUNT}-{ISOP_ABBR}{ISOX}'
                              .format(RING_TYPE=ring_type,
                                      SERIES_DB_COUNT=usr_db+isop_info_dct['DB'],
                                      ISOP_ABBR=isop_abbr,
                                      ISOX=class_type,
                                      )
                              )

            # for TxA, TxB
            else:
                abbr_formatter = ('{ISOP_SERIES}-{ISOP_ABBR}{ISOX}{SERIES_DB_COUNT}-'
                                  '{C_COUNT}:{DB_COUNT}[{DB_COUNT}xDB,'
                                  '{OH_COUNT}xOH,{KETO_COUNT}xKETO]'
                                  '<CHO@C{CHO_COUNT},COOH@{COOH_COUNT}>{left}OAP:{OAP_COUNT},OCP:{OCP_COUNT}{right}'
                                  .format(ISOP_SERIES=usr_series,
                                          RING_TYPE=ring_type,
                                          SERIES_DB_COUNT=usr_db+isop_info_dct['DB'],
                                          ISOP_ABBR=isop_abbr,
                                          ISOX=class_type,
                                          C_COUNT=usr_c,
                                          DB_COUNT=usr_db+isop_info_dct['DB'],
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
                abbr_short = ('{ISOP_ABBR}{ISOX}{SERIES_DB_COUNT}'
                              .format(SERIES_DB_COUNT=usr_db+isop_info_dct['DB'],
                                      ISOP_ABBR=isop_abbr,
                                      ISOX=class_type,
                                      )
                              )
        else:
            print(usr_c, usr_db+isop_info_dct['DB'], self.db_info_dct['DB_count'], 'not listed')
            abbr_formatter = ''
            abbr_short = ''

        return abbr_formatter, abbr_short

    def get_info_dct(self, isop_smi, isop_info_dct, usr_db, usr_pre_db_count, usr_pre_smi, reverse=False):

        isop_series = self.get_isop_series(isop_info_dct, usr_pre_smi, usr_pre_db_count, reverse=reverse)

        c_count = isop_smi.count('C')

        isop_checker, abbr_short = self.get_isop_checker(isop_info_dct, usr_c=c_count,
                                                         usr_db=usr_db, usr_series=isop_series)
        abbr_gen = AbbrGenerator()
        isop_abbr, isop_typ_str = abbr_gen.decode(isop_checker)

        isop_json = ('{left}"C": {C_COUNT}, "DB": {DB_COUNT}, "CHO": {CHO_COUNT}, "EPOXY": 0, '
                     '"OAP": {OAP_COUNT}, "OCP": {OCP_COUNT}, "COOH": {COOH_COUNT}, '
                     '"KETO": {KETO_COUNT}, "OH": {OH_COUNT}, "OOH": {OOH_COUNT}, '
                     '"LINK_TYPE": "", "Prostane": "{PROSTANE_TYPE}"{right}'
                     .format(PROSTANE_TYPE=abbr_short,
                             C_COUNT=c_count,
                             DB_COUNT=usr_db,
                             OH_COUNT=isop_info_dct['OH'],
                             OOH_COUNT=isop_info_dct['OOH'],
                             KETO_COUNT=isop_info_dct['KETO'],
                             CHO_COUNT=isop_info_dct['CHO'],
                             COOH_COUNT=isop_info_dct['COOH'],
                             left='{',
                             OAP_COUNT=isop_info_dct['OAP'],
                             OCP_COUNT=isop_info_dct['OCP'],
                             right='}'
                             )
                     )

        isop_lpp_sum_dct = {'SMILES': isop_smi, 'C_NUM': c_count, 'DB': usr_db,
                            'OAP': isop_info_dct['OAP'], 'OCP': isop_info_dct['OCP'],
                            'OH': isop_info_dct['OH'], 'KETO': isop_info_dct['KETO'], 'OOH': isop_info_dct['OOH'],
                            'CHO': isop_info_dct['CHO'], 'COOH': isop_info_dct['COOH'], 'EPOXY': 0,
                            'MOD_NUM': 1,
                            'FULL_SMILES': isop_smi,
                            'FA_CHECKER': isop_checker,
                            'FA_ABBR': isop_abbr, 'FA_TYPE': isop_typ_str, 'FA_JSON': isop_json,
                            'FRAG_SMILES': '[""]'}

        return isop_abbr, isop_lpp_sum_dct

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

            # print(_ring_info_dct['ORIGIN_SMILES'])
            # print(usr_db_main_smi)
            # _pattern_rgx = re.compile(_ring_info_dct['ORIGIN_SMILES'])

            _db_required_num = _ring_info_dct['DB_NUM_REQUIRED']
            _remain_db_num = usr_db_num - _db_required_num
            _pre_ring_db_count_lst = range(0, _remain_db_num + 1)
            _post_ring_db_count_lst = sorted(_pre_ring_db_count_lst, reverse=True)
            _pre_post_db_count_lst = zip(_pre_ring_db_count_lst, _post_ring_db_count_lst)
            # print(_pre_post_db_count_lst)
            for _db_seg in _pre_post_db_count_lst:
                # print('_db_seg', _db_seg)

                _rebuild_smi = 'C/C=C\\' * _db_seg[0] + _ring_info_dct['ORIGIN_SMILES'] + 'C/C=C\\' * _db_seg[1]
                # print('rebuild smi', _rebuild_smi)
                if _rebuild_smi == usr_db_main_smi:
                    # print('rebuild check pass!')
                    if _ring_info_dct['OCP'] == 0 and _ring_info_dct['OAP'] == 1:
                        _ring_lpp_smi = (usr_db_pre_smi + 'C/C=C\\' * _db_seg[0] +
                                         _ring_info_dct['MAIN_SMILES'] + 'C/C=C\\' * _db_seg[1] + usr_db_post_smi)

                        _isop_abbr, _isop_lpp_sum_dct = self.get_info_dct(_ring_lpp_smi, _ring_info_dct,
                                                                          _remain_db_num, _db_seg[0],
                                                                          usr_db_pre_smi, reverse=False)
                        # print('_ring_lpp_smi_OAP', _ring_lpp_smi)
                        usr_lpp_dct[_isop_abbr] = _isop_lpp_sum_dct
                    elif _ring_info_dct['OCP'] == 1 and _ring_info_dct['OAP'] == 0:
                        _ring_lpp_smi = (usr_db_pre_smi + 'C/C=C\\' * _db_seg[0] +
                                         _ring_info_dct['MAIN_SMILES'] + usr_db_post_smi)
                        _isop_abbr, _isop_lpp_sum_dct = self.get_info_dct(_ring_lpp_smi, _ring_info_dct,
                                                                          _remain_db_num, _db_seg[0],
                                                                          usr_db_pre_smi, reverse=False)
                        usr_lpp_dct[_isop_abbr] = _isop_lpp_sum_dct
                        # print('_ring_lpp_smi_OAP', _ring_lpp_smi)

                    _rev_smi = _ring_info_dct['REVERSE_SMILES']
                    if isinstance(_rev_smi, str):
                        # blank 'REVERSE_SMILES' gives np.nan which is float
                        if _rev_smi not in no_rev_lst:
                            if _ring_info_dct['OCP'] == 0 and _ring_info_dct['OAP'] == 1:
                                _ring_lpp_smi = (usr_db_pre_smi + 'C/C=C\\' * _db_seg[0] +
                                                 _ring_info_dct['REVERSE_SMILES']
                                                 + 'C/C=C\\' * _db_seg[1] + usr_db_post_smi)

                                _isop_abbr, _isop_lpp_sum_dct = self.get_info_dct(_ring_lpp_smi, _ring_info_dct,
                                                                                  _remain_db_num, _db_seg[0],
                                                                                  usr_db_pre_smi, reverse=True)
                                usr_lpp_dct[_isop_abbr] = _isop_lpp_sum_dct
                            elif _ring_info_dct['OCP'] == 1 and _ring_info_dct['OAP'] == 0:
                                _ring_lpp_smi = (usr_db_pre_smi + 'C/C=C\\' * _db_seg[0] +
                                                 _ring_info_dct['REVERSE_SMILES'] + usr_db_post_smi)
                                _isop_abbr, _isop_lpp_sum_dct = self.get_info_dct(_ring_lpp_smi, _ring_info_dct,
                                                                                  _remain_db_num, _db_seg[0],
                                                                                  usr_db_pre_smi, reverse=True)
                                usr_lpp_dct[_isop_abbr] = _isop_lpp_sum_dct
                        else:
                            pass
                    else:
                        pass

        # ox_isop_lpp_lst = [usr_db_pre_smi + _lpp_smi + usr_db_post_smi for _lpp_smi in usr_lpp_lst]
        # print(usr_lpp_dct.keys())

        return usr_lpp_dct


if __name__ == '__main__':

    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem

    isop_cfg = r'IsoP_ModConfig.csv'
    isop_abbr_cfg = r'IsoP_AbbrConfig.csv'

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

    got_lpp_dct = ox_isop.get_isop_lpp()

    lpp_mol_lst = []

    sdf_writer = Chem.SDWriter(save_sdf)

    '''
    'FULL_SMILES': isop_smi,
                            'FA_CHECKER': isop_checker,
                            'FA_ABBR': isop_abbr, 'FA_TYPE': isop_typ_str, 'FA_JSON': isop_json,
                            'FRAG_SMILES': '[""]'
    '''

    for _isop_lpp_abbr in got_lpp_dct.keys():
        _isop_lpp_info_dct = got_lpp_dct[_isop_lpp_abbr]
        _isop_lpp_smi = _isop_lpp_info_dct['FULL_SMILES']
        try:
            _mol = Chem.MolFromSmiles(_isop_lpp_smi)

            AllChem.Compute2DCoords(_mol)
            _mol.SetProp('_Name', _isop_lpp_info_dct['FA_ABBR'])
            _mol.SetProp('Description', _isop_lpp_info_dct['FA_ABBR'])
            _mol.SetProp('LM_ID', _isop_lpp_info_dct['FA_ABBR'])
            _mol.SetProp('ID', _isop_lpp_info_dct['FA_ABBR'])
            _mol.SetProp('FA_JSON', _isop_lpp_info_dct['FA_JSON'])
            _mol.SetProp('FA_ABBR', _isop_lpp_info_dct['FA_ABBR'])
            _mol.SetProp('SMILES', _isop_lpp_info_dct['FULL_SMILES'])
            sdf_writer.write(_mol)
            lpp_mol_lst.append(_mol)
        except:
            error_rgx = re.compile(r'//')
            error_checker = error_rgx.sub(r'', _isop_lpp_smi)
            # print(error_checker)
            try:
                _mol = Chem.MolFromSmiles(error_checker)
                AllChem.Compute2DCoords(_mol)
                lpp_mol_lst.append(_mol)
                print('corrected and pass')
            except:
                print('not working for --> ', _isop_lpp_smi)

                # pass

    # print(len(lpp_mol_lst))
    img = Draw.MolsToGridImage(lpp_mol_lst, molsPerRow=6, subImgSize=(200, 200))
    # img.show()
    img.save(save_img)
    print('saved!')
