# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import copy
import re
import pandas as pd

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
                          'PC_d9':
                              r'[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])(',
                          'd9-oxPC':
                              r'[O-]P(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])(OC[C@]([H])(',
                          'PC': r'[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(',
                          'PE': r'OP(OCCN)(OC[C@]([H])(',
                          'PG': r'OP(OCC(O)CO)(OC[C@]([H])(',
                          'PS': r'OP(OCC(C(O)=O)N)(OC[C@]([H])('}
        self.pl_hg_lst = ['PA', 'PC', 'PE', 'PG', 'PS', 'PC_d9', 'd9-oxPC']

        self.pl_end_smiles = r')=O'

        self.pl_hg_elem_dct = {}
        self.pl_hg_elem_dct['PC'] = {'hg': 'C5H14NO4P', 'hg_part': 'C3H9N'}
        self.pl_hg_elem_dct['d9-oxPC'] = {'hg': 'C5H5D9NO4P', 'hg_part': 'C3D9N'}

        scores_df = pd.read_csv('ion_scores.csv', index_col=0)
        # print scores_df.head()
        self.scores_dct = scores_df.to_dict()['ion_i']
        # print self.scores_dct

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
                    sn1_formula = s_obj.smiles2formula(fa_lst[1])
                    sn2_formula = s_obj.smiles2formula(fa_lst[0])
                    nl_formula_dct['sn1'] = sn1_formula
                    nl_formula_dct['sn2'] = sn2_formula
                    print sn1_formula, sn2_formula
                    print 'Get FA list', fa_lst
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

            theo_ion_mz_dct = {}

            for _ion in theo_ion_formula_dct.keys():

                _ion_formula = formula_obj.merge_dct(theo_ion_formula_dct[_ion], self.charge_dct[_chg])
                _ion_name_lst = ['[', _ion, _chg[2:]]
                _ion_name = ''.join(_ion_name_lst)
                theo_ion_mz_dct[_ion_name] = (mzcalc.get_mass(_ion_formula), self.scores_dct[_ion_name])

            print theo_ion_mz_dct.keys()
            _chg_info_dct['frag_info'] = theo_ion_mz_dct

            result_dct[_chg] = _chg_info_dct

        print result_dct

        sum_info_dct = {_name: result_dct}

        return sum_info_dct

    def frag2msp(self, sum_info_dct, outputname):

        output = open(outputname, mode='w')
        for _name in sum_info_dct.keys():
            info_dct = sum_info_dct[_name]
            for _chg in info_dct.keys():
                _chg_dct = info_dct[_chg]
                _frags_dct = _chg_dct['frag_info']
                _name_info = 'Name: ' + _name + ' ' + _name
                _id_info = 'LM_ID: ' + _name
                _pr_type = 'Precursor_type: ' + _chg
                _comment = 'Comment:' + str(_chg_dct['pr_mz'])
                _formula = _chg_dct['pr_formula']
                elem_count_checker = re.compile(r'([A-Z][a-z]|[A-Z])')
                elem_lst = elem_count_checker.split(_formula)
                # e.g. ['', 'C', '5', 'H', '14', 'N', '', 'O', '4', 'P', '', 'Na', '']
                for _i in range(len(elem_lst)):
                    if elem_lst[_i] == '1':
                        elem_lst[_i] = ''
                _p_formula = ''.join(elem_lst)
                _formula = 'Formula: ' + _p_formula
                _num_peak = 'Num Peaks: ' + str(len(_frags_dct.keys()))

                ion_lst = [_name_info, _id_info, _pr_type, _comment, _formula, _num_peak]
                for _ion in _frags_dct.keys():
                    _ion_name = '"'+ _ion + '"'
                    _ion_info_lst = [str(_frags_dct[_ion][0]), str(_frags_dct[_ion][1]), _ion_name]
                    # _ion_info_lst = [str(_frags_dct[_ion][0]), str(_frags_dct[_ion][1])]
                    _ion_info_txt = ' '.join(_ion_info_lst)
                    ion_lst.append(_ion_info_txt)
                ion_lst.append('\n')
                output.writelines('\n'.join(ion_lst))

        output.close()
        print 'msp Generated'


s = r'[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])(OC(CCCCCCC=O)=O)COC(CCCCCCCCCCCCCCC)=O)=O'

frag_obj = TheoFrag()

# dct = frag_obj.smiles2frag(s, 'test_LPP', chargelist=['[M+H]+', '[M+Na]+', '[M-H2O+H]+', '[M-H2O+Na]+'], plclass='PC')
dct = frag_obj.smiles2frag(s, 'test_LPP', chargelist=['[M-H]-', '[M+FA-H]-'], plclass='PC')

output = 'test.msp'

frag_obj.frag2msp(dct, output)
