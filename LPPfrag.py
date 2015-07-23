# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import copy
import re

from lpplibs.ExactMassCalc import Elem2Mass
from lpplibs.SMILESparser import SMILESparser
from lpplibs.FormulaCalc import FormulaCalc


class TheoFrag(object):

    def __init__(self):
        # define the charge type and elemental composition
        self.charge_dct = {'[M+H]+': {'H': 1}, '[M+Na]+': {'Na': 1}, '[M+K]+': {'K': 1},
                           '[M-H2O+H]+': {'H': -1, 'O': -1},
                           '[M-H2O+Na]+': {'H': -2, 'O': -1, 'Na': 1},
                           '[M-H2O+K]+': {'H': -2, 'O': -1, 'K': 1},
                           '[M-H]-': {'H': -1},
                           '[M+FA-H]-': {'C': 1, 'H': 1, 'O': 2},
                           '[M+HCOO]-': {'C': 1, 'H': 1, 'O': 2},
                           '[M-H2O-H]-': {'H': -3, 'O': -1},
                           '[M-H2O+FA-H]-': {'C': 1, 'H': -1, 'O': 1},
                           '[M-H2O+HCOO]-': {'C': 1, 'H': -1, 'O': 1}}
        self.charge_lst = ['[M+H]+', '[M+Na]+', '[M+K]+', '[M-H2O+H]+', '[M-H2O+Na]+', '[M-H2O+K]+',
                           '[M-H]-', '[M+FA-H]-', '[M+HCOO]-', '[M-H2O-H]-', '[M-H2O+FA-H]-', '[M-H2O+HCOO]-']

        self.pl_hg_dct = {'PA': r'OP(O)(OC[C@]([H])(',
                          'PC_d9': r'[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])(',
                          'PC': r'[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(',
                          'PE': r'OP(OCC[N])(OC[C@]([H])(',
                          'PG': r'OP(OCC(O)CO)(OC[C@]([H])(',
                          'PS': r'OP(OCC(C(O)=O)N)(OC[C@]([H])('}
        self.pl_hg_lst = ['PA', 'PC', 'PE', 'PG', 'PS', 'PC_d9']

        self.pl_end_smiles = r')=O'

        self.pl_hg_elem_dct = {}
        self.pl_hg_elem_dct['PC'] = {'hg': 'C5H14NO4P', 'hg_part': 'C3H9N'}

    def smiles2frag(self, usr_smiles, description, plclass=None, chargelist=None):

        """

        :param usr_smiles: str, The SMILES code to process
        :param description: str, The description of the SMILES code compound
        :param chargelist: list, list of charge status use the key in self.charge_dct e.g. ['[M+H]+', '[M+Na]+']
        :return: dict, result_dct
        """
        if not plclass:
            plclass = 'PC'
        if not chargelist:
            chargelist = []

        # load parameters
        _smiles = usr_smiles
        _name = description

        result_dct = {}  # e.g. result_dct['[M+H]+'] = {'pr_mz': 650, 'mz_df': mz_df}

        if plclass in self.pl_hg_lst:
            pl_hg_smiles = self.pl_hg_dct[plclass]
        else:
            result_dct['error'] = 'Lipid class error: ' + plclass
            return result_dct

        if chargelist == [] or chargelist == '':
            chargelist = ['[M+H]+', '[M+Na]+']
        else:
            for _chg in chargelist:
                if _chg in self.charge_lst:
                    pass
                else:
                    result_dct['error'] = 'Charge error: ' + _chg
                    return result_dct

        mzcalc = Elem2Mass()
        formula_obj = FormulaCalc()
        # start processing

        # Get neutral formula
        s_obj = SMILESparser()

        lpp_elem_str = s_obj.smiles2formula(_smiles)

        lpp_frag_dct = {}
        lpp_frag_dct['pr_formula'] = lpp_elem_str
        _lpp_formula_txt = mzcalc.get_elem(lpp_elem_str)
        _lpp_mass_f = mzcalc.get_mass(_lpp_formula_txt)
        lpp_frag_dct['pr_mass'] = _lpp_mass_f

        # process fragmentation
        nl_formula_dct = {}
        hg_elem_lst = self.pl_hg_elem_dct[plclass]
        for _hg in hg_elem_lst.keys():
            nl_formula_dct[_hg] = hg_elem_lst[_hg]

        _hg_length_i = len(pl_hg_smiles)

        if _smiles[:_hg_length_i] == pl_hg_smiles:
            _smiles_nohg = _smiles[_hg_length_i:]
            # Make re for connection between two FA
            fa_checker = re.compile(r'[)]COC[(]')

            if _smiles_nohg[-3:] == self.pl_end_smiles:
                # print _smiles_nohg[:-3]
                _smiles_sn1_sn2 = _smiles_nohg[:-3]
                fa_lst = fa_checker.split(_smiles_sn1_sn2)
                if len(fa_lst) == 2:
                    fa_lst[1] = 'OC(' + fa_lst[1]
                    sn1_formula = s_obj.smiles2formula(fa_lst[0])
                    sn2_formula = s_obj.smiles2formula(fa_lst[1])
                    nl_formula_dct['sn1'] = sn1_formula
                    nl_formula_dct['sn2'] = sn2_formula
                    # print sn1_formula, sn2_formula
                    # print 'Get FA list', fa_lst
                else:
                    result_dct['error'] = 'Can NOT process SMILES: ' + _smiles
            else:
                result_dct['error'] = 'Can NOT process SMILES: ' + _smiles
        else:
            result_dct['error'] = 'Can NOT process SMILES: ' + _smiles

        print nl_formula_dct

        theo_ion_formula_dct = copy.deepcopy(nl_formula_dct)

        for _nl in nl_formula_dct.keys():
            _ion_type = 'M-' + _nl
            theo_ion_formula_dct[_ion_type] = formula_obj.substract(lpp_elem_str, nl_formula_dct[_nl])

        # print theo_ion_formula_dct

        # print sn1_formula, sn2_formula
        # calc m/z for each charge type
        for _chg in chargelist:
            _chg_info_dct = {}
            _lpp_charged_str = s_obj.smiles2formula(_smiles, charge=_chg)
            _chg_info_dct['pr_formula'] = _lpp_charged_str

            _chg_formula_txt = mzcalc.get_elem(_lpp_charged_str)
            _chg_mz_f = mzcalc.get_mass(_chg_formula_txt)
            _chg_info_dct['pr_mz'] = _chg_mz_f

            result_dct[_chg] = _chg_info_dct

        print result_dct


s = r'[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(CCCCCCC=O)=O)COC(CCCCCCCCCCCCCCC)=O)=O'

frag_obj = TheoFrag()

frag_obj.smiles2frag(s, 'POPC', chargelist=['[M+H]+', '[M+Na]+', '[M+K]+'], plclass='PC')
