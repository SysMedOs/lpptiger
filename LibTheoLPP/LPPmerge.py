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
# import re

from rdkit import Chem
from rdkit.Chem import Draw


class LPPsmi:
    def __init__(self):
        print("Start to Merge back LPP-->")

    @staticmethod
    def smiles2abbr(usr_sn1, usr_sn2):
        print(usr_sn1)
        print(usr_sn2)

    @staticmethod
    def smiles2img(usr_smiles):
        _mol = Chem.MolFromSmiles(usr_smiles)
        _img = Draw.MolToImage(_mol, size=(900, 900))
        _img.show()
        return _mol

    @staticmethod
    def smiles2mol(usr_smiles):
        _mol = Chem.MolFromSmiles(usr_smiles)
        return _mol


# construct a decorator
def lpp_merge(theolpp_cls):
    """

    :param theolpp_cls:
    :return:
    """

    def _lpp_merge(px_lpp):
        def __lpp_merge(usr_hg, sn1=None, sn2=None):

            lpp_str = px_lpp(usr_hg, sn1, sn2)
            # print(lpp_str)
            # theolpp_cls.smiles2abbr(usr_sn1, usr_sn2)
            # theolpp_cls.smiles2img(lpp_str)

            return lpp_str
        return __lpp_merge
    return _lpp_merge


@lpp_merge(LPPsmi)
def pl_lpp(usr_hg, sn1=None, sn2=None):
    pl_hg_dct = {'PA': r'OP(O)(OCC(',
                 'PC': r'[O-]P(OCC[N+](C)(C)C)(OCC(',
                 'PC-CH3': r'[O]P(OCC[N](C)C)(OCC(',
                 'PE': r'OP(OCCN)(OCC(',
                 'PG': r'OP(OCC(O)CO)(OCC(',
                 'PS': r'OP(OCC(C(O)=O)N)(OCC(',
                 'PI': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O))(OCC(',
                 'PIP': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC(',
                 'PI4P': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(OCC('}

    if usr_hg.upper() in pl_hg_dct.keys():
        pl_hg = pl_hg_dct[usr_hg.upper()]
        gly_part = r')C'
        # sn1 FA 16:0
        # sn1 = r'OC(CCCCCCCCCCCCCCC)=O'
        # sn2 FA 18:2 (9Z, 12Z)
        # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
        pl_end = r')=O'
        # the following order is important!!
        pl_str = ''.join([pl_hg, sn2, gly_part, sn1, pl_end])
        return pl_str


@lpp_merge(LPPsmi)
def pl_hg_lpp(usr_hg, sn1=None, sn2=None):
    pl_hg_dct = {'PA': r'OP(O)(O)=O',
                 'PC': r'[O-]P(OCC[N+](C)(C)C)(O)=O',
                 'PE': r'OP(OCCN)(O)=O',
                 'PG': r'OP(OCC(O)CO)(O)=O',
                 'PS': r'OP(OCC(C(O)=O)N)(O)=O',
                 'PI': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O))(O)=O',
                 'PIP': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(O)=O',
                 'PI4P': r'OP(O[C@H]1[C@H](O)([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O))(O)=O'}

    gly_pre_part = r'/C=C/(',

    if usr_hg.upper() in pl_hg_dct.keys():
        pl_hg = pl_hg_dct[usr_hg.upper()]
        gly_part = r')C'
        # sn1 FA 16:0
        # sn1 = r'OC(CCCCCCCCCCCCCCC)=O'
        # sn2 FA 18:2 (9Z, 12Z)
        # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
        pl_end = r')=O'
        # the following order is important!!
        pl_str = ''.join([pl_hg, sn2, gly_part, sn1, pl_end])
        return pl_str
