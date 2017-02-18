# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re


class AbbrDecoder(object):

    def __init__(self):
        pass

    @staticmethod
    def get_main_parts(usr_abbr):

        """
        destinguish PL abbr.
        :param usr_abbr: (str) user input
        :return:
        """

        pl_full_rgx = re.compile(r'(P[ACEGIS])(\()(.{1,})(\))')
        pl_rest_rgx = re.compile(r'(.*)(/)(.*)')

        _pl_checker = re.match(pl_full_rgx, usr_abbr)
        try:
            if _pl_checker:
                pre_pl_lst = _pl_checker.groups()
                pl_hg = pre_pl_lst[0]

                pl_rest = pre_pl_lst[2]

                _pl_rest_check = re.match(pl_rest_rgx, pl_rest)

                if _pl_rest_check:
                    rest_lst = _pl_rest_check.groups()
                    pl_sn1 = rest_lst[0]
                    pl_sn2 = rest_lst[2]
                    return pl_hg, pl_sn1, pl_sn2

        except:
            _msg = 'PL abbr %s is not correct ' % usr_abbr
            return _msg, 'error'

    def pl_class(self):
        def calc_elem(func):
            _x = ''.join(func)
            z = [_x]
            print 'z', z
            return z

    @pl_class
    def get_pl_elem(self, usr_hg):
        _a = list(usr_hg)
        return _a




# t1 = 'PA(16:0/18:1)'
# t2 = 'PAx(0:0/18:1)'
# t3 = 'PG(16:0/18:2)'
# t3 = 'PC(18:0/20:4)'
# l = [t1, t2, t3]
# a = AbbrDecoder()
# for i in l:
#     x = a.get_main_parts(i)
#     print x
#     b = x[0]
#     print b
#     m = a.get_pl_elem(b)
# print 'fin!'



