# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re


class AbbrGenerator(object):

    def __init__(self):
        pass

    def decode(self, usr_code):

        # prepare output
        lpp_abbr_str = ''
        lpp_typ_str = ''

        # usr_code = 'P-18:0[0xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'

        lpp_code_rgx = re.compile(r'(?P<fa>\d{1,2}[:]\d|[OPop]-\d{1,2}[:]\d)([\[])(?P<mod>.*)([\]][<])'
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

# x = 'P-18:1[1xDB,0xOH,1xKETO]<CHO@C0,COOH@C0>{OAP:1,OCP:0}'
# x = 'P-18:0[0xDB,0xOH,0xKETO]<CHO@C0,COOH@C0>{OAP:0,OCP:0}'
# x = '9:0[0xDB,1xOH,0xKETO]<CHO@C9,COOH@C0>{OAP:0,OCP:1}'
#
# a = AbbrGenerator()
# m, n = a.decode(x)
#
# print(m, n)
