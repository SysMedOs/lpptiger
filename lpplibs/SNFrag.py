# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
import re

from rdkit import Chem

from DBoxTheo import fa_link_filter


class SNFrag(object):

    def __abs__(self):
        pass

    def get_water_loss(self, usr_sn_smiles):

        _sn_link_dct = fa_link_filter(usr_sn_smiles)

        print(_sn_link_dct)

        _sn_smiles = _sn_link_dct['FULL_smiles']
        _sn_pre_str = _sn_link_dct['PRE_str']
        _sn_post_str = _sn_link_dct['POST_str']
        _sn_pre_rgx_str = _sn_link_dct['PRE_rgx']
        _sn_post_rgx_str = _sn_link_dct['POST_rgx']
        _sn_link_str = _sn_link_dct['LINK_type']

        # create output list
        sn_smi_lst = [_sn_smiles]

        sn_rgx = re.compile(r'(%s)(.*)(%s)' % (_sn_pre_rgx_str, _sn_post_rgx_str))
        sn_hydro_rgx = re.compile(r'([(]O[)])')
        # before or after OH maybe C=C
        sn_pre_db_rgx = re.compile(r'(.*)([/\\]?C[=]C[/\\]?C)')
        sn_post_db_rgx = re.compile(r'([/\\]?C[=]C[/\\]?)(.*)')

        sn_checker = re.match(sn_rgx, _sn_smiles)
        if sn_checker:
            pre_sn_lst = sn_checker.groups()
            print('pre_sn_lst', pre_sn_lst)
            sn_main_str = pre_sn_lst[1]
            sn_hydro_checker = re.split(sn_hydro_rgx, sn_main_str)
            if sn_hydro_checker:
                pre_hydro_lst = sn_hydro_checker
                print('pre_hydro_lst', pre_hydro_lst)
                if len(pre_hydro_lst) >= 3:
                    _sn_pre_hydro_smi = pre_hydro_lst[0]
                    _sn_post_hydro_smi = pre_hydro_lst[2]

                    print('_sn_pre_hydro_smi', _sn_pre_hydro_smi)
                    print('_sn_post_hydro_smi', _sn_post_hydro_smi)

                    # sn_pre_db_checker = re.match(sn_pre_db_rgx, _sn_pre_hydro_smi)
                    # sn_post_db_checker = re.match(sn_post_db_rgx, _sn_post_hydro_smi)
                    #
                    # if sn_pre_db_checker and not sn_post_db_checker:
                    #     sn_post_db_lst = sn_post_db_checker.groups()

        else:
            print('no Hydroxyl groups')

usr_smi = r'OC(CCCCC(O)/C=C/C(O)/C=C/C(O)/C=C/CC(O)=O)=O'

a = SNFrag()

a.get_water_loss(usr_smi)
