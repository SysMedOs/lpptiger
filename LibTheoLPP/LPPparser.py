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

import pandas as pd


class PLParser(object):

    def __init__(self):
        pass

    @staticmethod
    def check_lyso(usr_pl_abbr):
        lyso_rgx = re.compile(r'[Li][Yy][Ss][Oo]')
        lyso_groups_rgx = re.compile(r'([Li][Yy][Ss][Oo])(P[ACEGIS]|PI\d?P)'
                                     r'([(])(\d\d\s?[:]\s?\d|0[:]0)(\s?[/_;]\s?)(\d\d\s?[:]\s?\d)(\s?[)])')
        lyso_checker = re.match(lyso_rgx, usr_pl_abbr)
        if lyso_checker:
            pass

    @staticmethod
    def get_composition(usr_pl_abbr):
        pl_groups_rgx = re.compile(r'(?P<hg>P[ACEGIS]|PI\d?P)(\s?[(]\s?)'
                                   r'(?P<sn1>\d\d\s?[:]\s?\d[n\-\d]{0,4}|[OP]\s?-\s?\d\d\s?[:]\s?\d[n\-\d]{0,4})'
                                   r'(\s?[/_;]\s?)'
                                   r'(?P<sn2>\d\d\s?[:]\s?\d[n\-\d]{0,4}|[OP]\s?-\s?\d\d\s?[:]\s?\d[n\-\d]{0,4})'
                                   r'(\s?[)])')
        space_rgx = re.compile('\s+')

        pl_groups_checker = re.match(pl_groups_rgx, usr_pl_abbr)
        if pl_groups_checker:
            # remove all spaces in the elements by re and lambda function, get only info from hg, sn1 and sn2
            pre_pl_elem_lst = map(lambda x: re.sub(space_rgx, '', x), pl_groups_checker.group('hg', 'sn1', 'sn2'))
        else:
            pre_pl_elem_lst = []

        if len(pre_pl_elem_lst) > 0:

            pl_elem_lst = [pre_pl_elem_lst[0]]
            pl_info_dct = {'HG': pre_pl_elem_lst[0]}
            
            fa_omega_rgx = re.compile(r'(?P<linker>[OP]\s?-\s?|\s?)(?P<carbon>\d\d)(\s?[:]\s?)'
                                      r'(?P<db>\d)(?P<omega>n-\d{1,2}|\s?)')
            omega_rgx = re.compile(r'(n-)(?P<position>\d{1,2})')
            
            sn1_abbr = pre_pl_elem_lst[1]
            
            sn1_checker = re.match(fa_omega_rgx, sn1_abbr)
            if sn1_checker:
                sn1_link = sn1_checker.group('linker')
                sn1_c_num = sn1_checker.group('carbon')
                sn1_db_num = sn1_checker.group('db')
                sn1_omega_abbr = sn1_checker.group('omega')
                if len(sn1_omega_abbr) > 0:
                    sn1_omega_checker = re.match(omega_rgx, sn1_omega_abbr)
                    if sn1_omega_checker:
                        sn1_omega_type = sn1_omega_checker.group('position')
                    else:
                        sn1_omega_type = '0'
                else:
                    sn1_omega_type = '0'
            else:
                sn1_link = ''
                sn1_c_num = ''
                sn1_db_num = ''
                sn1_omega_type = '0'

            sn2_abbr = pre_pl_elem_lst[2]
            sn2_checker = re.match(fa_omega_rgx, sn2_abbr)
            if sn2_checker:
                sn2_link = sn2_checker.group('linker')
                sn2_c_num = sn2_checker.group('carbon')
                sn2_db_num = sn2_checker.group('db')
                sn2_omega_abbr = sn2_checker.group('omega')
                if len(sn2_omega_abbr) > 0:
                    sn2_omega_checker = re.match(omega_rgx, sn2_omega_abbr)
                    if sn2_omega_checker:
                        sn2_omega_type = sn2_omega_checker.group('position')
                    else:
                        sn2_omega_type = '0'
                else:
                    sn2_omega_type = '0'
            else:
                sn2_link = ''
                sn2_c_num = ''
                sn2_db_num = ''
                sn2_omega_type = '0'

            pl_elem_lst.append(''.join([sn1_link, sn1_c_num, ':', sn1_db_num]))
            pl_elem_lst.append(''.join([sn2_link, sn2_c_num, ':', sn2_db_num]))

            link_dct = {'O-': 'O', 'P-': 'P', '': 'A'}

            pl_info_dct['sn1_link'] = link_dct[sn1_link]
            pl_info_dct['sn2_link'] = link_dct[sn2_link]
            pl_info_dct['sn1_c_num'] = sn1_c_num
            pl_info_dct['sn2_c_num'] = sn2_c_num
            pl_info_dct['sn1_db_num'] = sn1_db_num
            pl_info_dct['sn2_db_num'] = sn2_db_num
            pl_info_dct['sn1_omega_type'] = sn1_omega_type
            pl_info_dct['sn2_omega_type'] = sn2_omega_type
            pl_info_dct['sn1_abbr'] = pl_elem_lst[1]
            pl_info_dct['sn2_abbr'] = pl_elem_lst[2]

        else:
            pl_elem_lst = []
            pl_info_dct = {}

        return pl_elem_lst, pl_info_dct
