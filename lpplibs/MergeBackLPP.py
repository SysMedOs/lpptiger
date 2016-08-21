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
        def __lpp_merge(usr_sn1, usr_sn2):

            lpp_str = px_lpp(usr_sn1, usr_sn2)
            print(lpp_str)
            theolpp_cls.smiles2abbr(usr_sn1, usr_sn2)
            theolpp_cls.smiles2img(lpp_str)

        return __lpp_merge
    return _lpp_merge


@lpp_merge(TheoLPP)
def pa_lpp(usr_sn1, usr_sn2):
    pa_hg = r'OP(O)(OCC([H])('
    gly_part = r')C'
    # sn1 FA 16:0
    # sn1 = r'OC(CCCCCCCCCCCCCCC)=O'
    # sn2 FA 18:2 (9Z, 12Z)
    # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
    pl_end = r')=O'
    # the following order is important!!
    pl_str = ''.join([pa_hg, usr_sn1, gly_part, usr_sn2, pl_end])
    return pl_str


@lpp_merge(TheoLPP)
def pc_lpp(usr_sn1, usr_sn2):
    pc_hg = r'[O-]P(OCC[N+](C)(C)C)(OCC([H])('
    # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
    gly_part = r')C'
    # sn2 = r'OC(CCCCCCCCCCCCCCC)=O'
    pl_end = r')=O'
    pl_str = ''.join([pc_hg, usr_sn1, gly_part, usr_sn2, pl_end])
    return pl_str


@lpp_merge(TheoLPP)
def pe_lpp(usr_sn1, usr_sn2):
    pe_hg = r'OP(OCCN)(OCC([H])('
    # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
    gly_part = r')C'
    # sn2 = r'OC(CCCCCCCCCCCCCCC)=O'
    pl_end = r')=O'
    pl_str = ''.join([pe_hg, usr_sn1, gly_part, usr_sn2, pl_end])
    return pl_str


@lpp_merge(TheoLPP)
def pg_lpp(usr_sn1, usr_sn2):
    pg_hg = r'OP(OCC(O)CO)(OCC([H])('
    # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
    gly_part = r')C'
    # sn2 = r'OC(CCCCCCCCCCCCCCC)=O'
    pl_end = r')=O'
    pl_str = ''.join([pg_hg, usr_sn1, gly_part, usr_sn2, pl_end])
    return pl_str


@lpp_merge(TheoLPP)
def ps_lpp(usr_sn1, usr_sn2):
    ps_hg = r'OP(OCC(C(O)=O)N)(OCC([H])('
    # sn2 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
    gly_part = r')C'
    # sn2 = r'OC(CCCCCCCCCCCCCCC)=O'
    pl_end = r')=O'
    pl_str = ''.join([ps_hg, usr_sn1, gly_part, usr_sn2, pl_end])
    return pl_str


sn1 = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
sn2 = r'OC(CCCCCCCCCCCCCCC)=O'

pa_lpp(sn1, sn2)
pc_lpp(sn1, sn2)
pe_lpp(sn1, sn2)
pg_lpp(sn1, sn2)
ps_lpp(sn1, sn2)
