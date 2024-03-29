# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2017  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Maria Fedorova
# LPPtiger is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
#
# For more info please contact:
#     SysMedOs_team: oxlpp@bbz.uni-leipzig.de
#     LPPtiger repository: https://bitbucket.org/SysMedOs/lpptiger
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#

from __future__ import print_function
import re
import json

import  pandas as pd


def fa_abbr_encode(fa_dct, mod_info_t_df, fa_ox_df, mode=0):

    fa_link_typ = fa_dct['DB_LINK_type']
    db_count = fa_dct['DB_count']
    db_start_count = fa_dct['DB_start']

    db_pre_part = fa_dct['DB_pre_part']
    db_post_part = fa_dct['DB_post_part']

    db_count_i_lst = range(1, db_count + 1)
    db_count_s_lst = [str(db_x) for db_x in db_count_i_lst]

    mod_smiles_lst = mod_info_t_df['SMILES'].tolist()
    mod_dct = {}
    for _idx, _mod_r in mod_info_t_df.iterrows():
        # after .t and reset_index, the mod name became column named 'index'
        mod_dct[_mod_r['SMILES']] = _idx

    fa_ox_df['FULL_SMILES'] = db_pre_part + fa_ox_df[db_count_s_lst].sum(axis=1) + db_post_part

    fa_ox_df['ABBR'] = ''
    left_dct_str = '{'
    right_dct_str = '}'

    for _idx, fa_r in fa_ox_df.iterrows():
        _db_counter = 0
        _oh_counter = 0
        _ooh_counter = 0
        _keto_counter = 0
        _epoxy_counter = 0
        _cho_counter = 0
        _cooh_counter = 0
        _oap_counter = 0
        _ocp_counter = 0
        _fa_smi = fa_r['FULL_SMILES']
        _c_counter = _fa_smi.count('C')
        for db_s in db_count_s_lst:
            db_smi = fa_r[db_s]
            print(db_smi)
            if db_smi in mod_smiles_lst:
                mod_idx = mod_dct[db_smi]
                _db_counter += mod_info_t_df.get_value(mod_idx, 'DB')
                _oh_counter += mod_info_t_df.get_value(mod_idx, 'OH')
                _ooh_counter += mod_info_t_df.get_value(mod_idx, 'OOH')
                _keto_counter += mod_info_t_df.get_value(mod_idx, 'KETO')
                _epoxy_counter += mod_info_t_df.get_value(mod_idx, 'EPOXY')
                _cho_counter += mod_info_t_df.get_value(mod_idx, 'CHO')
                _cooh_counter += mod_info_t_df.get_value(mod_idx, 'COOH')
                _oap_counter += mod_info_t_df.get_value(mod_idx, 'OAP')
                _ocp_counter += mod_info_t_df.get_value(mod_idx, 'OCP')
            elif db_smi == 'C/C=C/':
                _db_counter += 1
        print(_fa_smi)
        print(_db_counter, _oh_counter, _ooh_counter, _keto_counter)
        if mode == 0:
            if _ocp_counter == 0:
                _fa_abbr_str = ('{LINK_TYPE}{NUM_C}:{NUM_DB}'
                                '[{NUM_DB}xDB,{NUM_OH}xOH,{NUM_KETO}xKETO,{NUM_OOH}xOOH,{NUM_EPOXY}xEPOXY]'
                                .format(LINK_TYPE=fa_link_typ,
                                        NUM_C=_c_counter,
                                        NUM_DB=_db_counter,
                                        NUM_OH=_oh_counter,
                                        NUM_KETO=_keto_counter,
                                        NUM_OOH=_ooh_counter,
                                        NUM_EPOXY=_epoxy_counter
                                        )
                                )
                fa_ox_df.set_value(_idx, 'ABBR', _fa_abbr_str)
            elif _ocp_counter == 1:
                if _cho_counter == 1 and _cooh_counter == 0:
                    _fa_abbr_str = ('{LINK_TYPE}{NUM_C}:{NUM_DB}'
                                    '[{NUM_DB}xDB,{NUM_OH}xOH,{NUM_KETO}xKETO,{NUM_OOH}xOOH,{NUM_EPOXY}xEPOXY]'
                                    '<CHO@C{NUM_C}>'
                                    .format(LINK_TYPE=fa_link_typ,
                                            NUM_C=_c_counter,
                                            NUM_DB=_db_counter,
                                            NUM_OH=_oh_counter,
                                            NUM_KETO=_keto_counter,
                                            NUM_OOH=_ooh_counter,
                                            NUM_EPOXY=_epoxy_counter
                                            )
                                    )
                    fa_ox_df.set_value(_idx, 'ABBR', _fa_abbr_str)
                if _cho_counter == 0 and _cooh_counter == 1:
                    _fa_abbr_str = ('{LINK_TYPE}{NUM_C}:{NUM_DB}'
                                    '[{NUM_DB}xDB,{NUM_OH}xOH,{NUM_KETO}xKETO,{NUM_OOH}xOOH,{NUM_EPOXY}xEPOXY]'
                                    '<COOH@C{NUM_C}>'
                                    .format(LINK_TYPE=fa_link_typ,
                                            NUM_C=_c_counter,
                                            NUM_DB=_db_counter,
                                            NUM_OH=_oh_counter,
                                            NUM_KETO=_keto_counter,
                                            NUM_OOH=_ooh_counter,
                                            NUM_EPOXY=_epoxy_counter
                                            )
                                    )
                    fa_ox_df.set_value(_idx, 'ABBR', _fa_abbr_str)
            else:
                pass
        else:
            # todo zhixu.ni@uni-leipzig.de: add exact position of modification sites e.g. 12:1[1xDB{7},1xOH{9}]

            _fa_abbr_str = ('{LINK_TYPE}{NUM_C}:{NUM_DB}'
                            '[{NUM_DB}xDB,{NUM_OH}xOH,{NUM_KETO}xKETO,{NUM_OOH}xOOH,{NUM_EPOXY}xEPOXY]'
                            .format(LINK_TYPE=fa_link_typ,
                                    NUM_C=_c_counter,
                                    NUM_DB=_db_counter,
                                    NUM_OH=_oh_counter,
                                    NUM_KETO=_keto_counter,
                                    NUM_OOH=_ooh_counter,
                                    NUM_EPOXY=_epoxy_counter
                                    )
                            )
            fa_ox_df.set_value(_idx, 'ABBR', _fa_abbr_str)

    return fa_ox_df


class AbbrGenerator(object):

    def __init__(self):
        pass

    @staticmethod
    def decode(usr_code):

        # prepare output
        lpp_abbr_str = ''
        lpp_typ_str = ''

        # usr_code = 'P-18:0[0xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'

        # only normal FA and O- P-
        # lpp_code_rgx = re.compile(r'(?P<fa>\d{1,2}[:]\d|[OPop]-\d{1,2}[:]\d)([\[])(?P<mod>.*)([\]][<])'
        #                           r'(?P<end>.*)([>][{])(?P<typ>.*)([}])')

        # Also for IsoP
        lpp_code_rgx = re.compile(r'(?P<fa>\d{1,2}[:]\d|.*-\d{1,2}[:]\d)([\[])(?P<mod>.*)([\]][<])'
                                  r'(?P<end>.*)([>][{])(?P<typ>.*)([}])')

        mod_rgx = re.compile(r'[1-9]x\w{2,13}')  # filter out 0 mods, support to hydroperoxyl (12 letters)
        end_rgx = re.compile(r'CHO@C1|COOH@C1')  # filter out the correct end
        typ_rgx = re.compile(r'O[AC]P[:][1-9]')

        lpp_code_checker = re.match(lpp_code_rgx, usr_code)
        if lpp_code_checker:
            fa_info_dct = lpp_code_checker.groupdict()  # dict for 'fa', 'mod', 'end', 'typ'

            mod_lst = re.findall(mod_rgx, fa_info_dct['mod'])
            end_lst = re.findall(end_rgx, fa_info_dct['end'])
            typ_lst = re.findall(typ_rgx, fa_info_dct['typ'])

            # print(mod_lst)
            # print(end_lst)
            # print(typ_lst)

            # check if unmodified FA
            if len(mod_lst) + len(end_lst) == 0 and len(typ_lst) == 0:
                if fa_info_dct['fa'][-2:] == ':0':
                    lpp_abbr_str = fa_info_dct['fa']
                    lpp_typ_str = ''
            # check if OAP or OCP
            if len(mod_lst) + len(end_lst) > 0 and len(typ_lst) >= 1:

                if len(typ_lst) == 1:
                    lpp_typ_str = typ_lst[0][0:3]
                elif len(typ_lst) == 2:
                    lpp_typ_str = 'OCP'
                else:
                    lpp_typ_str = ''

                if len(mod_lst) > 0:
                    _mod_str = ','.join(mod_lst)
                    _mod_str = ''.join(['[', _mod_str, ']'])
                else:
                    _mod_str = ''
                if len(end_lst) == 1:
                    _end_str = ','.join(end_lst)
                    _end_str = ''.join(['<', _end_str[:-1], fa_info_dct['fa'].split(':')[0], '>'])
                else:
                    _end_str = ''

                lpp_abbr_str = ''.join([fa_info_dct['fa'], _mod_str, _end_str])

            # lpp_info_lst = (lpp_abbr_str, lpp_typ_str)
            return lpp_abbr_str, lpp_typ_str

    @staticmethod
    def from_json(usr_json_str):

        # prepare output
        lpp_abbr_str = ''
        lpp_typ_str = ''

        fa_dct = json.load(usr_json_str)
        if fa_dct['OAP'] == 1 and fa_dct['OCP'] == 0:
            lpp_typ_str = 'OAP'
        if fa_dct['OAP'] == 0 and fa_dct['OCP'] == 1:
            lpp_typ_str = 'OCP'
        if fa_dct['OAP'] == 0 and fa_dct['OCP'] == 0:
            lpp_typ_str = ''

        lpp_abbr_str = ''.join([fa_dct['C_num'], ':', fa_dct['DB_num']])
        if fa_dct['MOD_NUM'] > 0:
            lpp_mod_lst = []


# x = 'P-18:1[1xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'
# x = 'P-18:0[0xDB,0xOH,0xKETO]<CHO@C0,COOH@C0>{OAP:0,OCP:0}'
# x = '9:0[0xDB,1xOH,0xKETO]<CHO@C9,COOH@C0>{OAP:0,OCP:1}'
#
# a = AbbrGenerator()
# m, n = a.decode(x)
#
# print(m, n)
