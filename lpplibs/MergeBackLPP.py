# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de


from __future__ import print_function
# import re

from rdkit import Chem
from rdkit.Chem import Draw


class TheoLPP:
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


@lpp_merge(TheoLPP)
def pl_lpp(usr_hg, sn1=None, sn2=None):
    pl_hg_dct = {'PA': r'OP(O)(OCC(',
                 'PC': r'[O-]P(OCC[N+](C)(C)C)(OCC(',
                 'PE': r'OP(OCCN)(OCC(',
                 'PG': r'OP(OCC(O)CO)(OCC(',
                 'PS': r'OP(OCC(C(O)=O)N)(OCC(',
                 'PI': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[O])(OCC(',
                 'PIP': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O])(OCC(',
                 'PI4P': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O])(OCC('}

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


@lpp_merge(TheoLPP)
def pl_hg_lpp(usr_hg, sn1=None, sn2=None):
    pl_hg_dct = {'PA': r'OP(O)(O)=O',
                 'PC': r'[O-]P(OCC[N+](C)(C)C)(O)=O',
                 'PE': r'OP(OCCN)(O)=O',
                 'PG': r'OP(OCC(O)CO)(O)=O',
                 'PS': r'OP(OCC(C(O)=O)N)(O)=O',
                 'PI': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[O])(O)=O',
                 'PIP': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O])(O)=O',
                 'PI4P': r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O])(O)=O'}

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


# @lpp_merge(TheoLPP)
# def pa_lpp(usr_sn1, usr_sn2):
#     pa_hg = r'OP(O)(OCC([H])('
#     gly_part = r')C'
#     # sn1 FA 16:0
#     # sn1 = r'OC(CCCCCCCCCCCCCCC)=O'
#     # sn2 FA 18:2 (9Z, 12Z)
#     # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
#     pl_end = r')=O'
#     # the following order is important!!
#     pl_str = ''.join([pa_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def pc_lpp(usr_sn1, usr_sn2):
#     pc_hg = r'[O-]P(OCC[N+](C)(C)C)(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([pc_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def pe_lpp(usr_sn1, usr_sn2):
#     pe_hg = r'OP(OCCN)(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([pe_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def pg_lpp(usr_sn1, usr_sn2):
#     pg_hg = r'OP(OCC(O)CO)(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([pg_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def pi_lpp(usr_sn1, usr_sn2):
#     # Inositol, IUPAC (1R,2R,3S,4S,5R,6S)-cyclohexane-1,2,3,4,5,6-hexol
#     # O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[O]
#     pi_hg = r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[O])(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([pi_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def pi4p_lpp(usr_sn1, usr_sn2):
#     # Inositol-4-phospate
#     # O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O]
#     pi4p_hg = r'OP(O[C@H]1[C@@H]([C@@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H]1O)[O])(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([pi4p_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str
# 
# 
# @lpp_merge(TheoLPP)
# def ps_lpp(usr_sn1, usr_sn2):
#     ps_hg = r'OP(OCC(C(O)=O)N)(OCC([H])('
#     gly_part = r')C'
#     pl_end = r')=O'
#     pl_str = ''.join([ps_hg, usr_sn1, gly_part, usr_sn2, pl_end])
#     return pl_str


# sn1 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
# sn2 = r'OC(CCCCCCCCCCCCCCC)=O'
#
# # pa_lpp(sn1, sn2)
# # pc_lpp(sn1, sn2)
# # pe_lpp(sn1, sn2)
# # pg_lpp(sn1, sn2)
# # pi_lpp(sn1, sn2)
# # pi4p_lpp(sn1, sn2)
# # ps_lpp(sn1, sn2)
#
# x = pc_lpp(sn1, sn2)
# print('x', x)
