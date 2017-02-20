# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of TheoLPP.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de


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
                                   r'(?P<sn1>\d\d\s?[:]\s?\d|[OP]\s?-\s?\d\d\s?[:]\s?\d)(\s?[/_;]\s?)'
                                   r'(?P<sn2>\d\d\s?[:]\s?\d|[OP]\s?-\s?\d\d\s?[:]\s?\d)(\s?[)])')
        space_rgx = re.compile('\s+')

        pl_groups_checker = re.match(pl_groups_rgx, usr_pl_abbr)
        if pl_groups_checker:
            # remove all spaces in the elements by re and lambda function, get only info from hg, sn1 and sn2
            pl_elem_lst = map(lambda x: re.sub(space_rgx, '', x), pl_groups_checker.group('hg', 'sn1', 'sn2'))
        else:
            pl_elem_lst = []

        return pl_elem_lst

