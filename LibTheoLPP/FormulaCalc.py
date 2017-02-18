# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
import copy
from ExactMassCalc import Elem2Mass


class FormulaCalc(object):

    def formula2dct(self, formula):

        elem_whitelist = ['H', 'D', 'C', 'O', 'P', 'N', 'Na', 'K', 'S', 'Cl', 'Ca', 'Fe', 'Cu']

        elem_count_checker = re.compile(r'([A-Z][a-z]|[A-Z])')
        elem_lst = elem_count_checker.split(formula)
        # e.g. ['', 'C', '5', 'H', '14', 'N', '', 'O', '4', 'P', '', 'Na', '']
        if elem_lst[0] == '':
            elem_lst.pop(0)
            # e.g. ['C', '5', 'H', '14', 'N', '', 'O', '4', 'P', '', 'Na', '']

        elem_dct = {}
        idx_lst = range(len(elem_lst))

        for idx in idx_lst:
            _atom = elem_lst[idx]
            try:
                _count = elem_lst[idx + 1]
                if _atom in elem_whitelist:
                    if _count == '':
                        elem_dct[_atom] = 1
                    else:
                        elem_dct[_atom] = int(_count)
                else:
                    pass
            except (IndexError, ValueError):
                pass

        return elem_dct

    def merge(self, _usr_formula_lst=None):

        formula_sum_str = ''
        formula_sum_dct = {}

        if not _usr_formula_lst:
            formula_sum_str = 'Error: No formula list'
            return formula_sum_str

        if isinstance(_usr_formula_lst, list) and len(_usr_formula_lst) == 1:
            formula_sum_str = 'Error: need more than one formula'
            return formula_sum_str

        elif isinstance(_usr_formula_lst, list) and len(_usr_formula_lst) > 1:
            masscalc = Elem2Mass()
            parsed_formula_lst = []
            for _formula in _usr_formula_lst:
                _parsed_formula = masscalc.get_elem(_formula)
                parsed_formula_lst.append(_parsed_formula)

            for _p_formula in parsed_formula_lst:
                _p_dct = self.formula2dct(_p_formula)
                for _elem in _p_dct.keys():
                    if _elem in formula_sum_dct.keys():
                        formula_sum_dct[_elem] += _p_dct[_elem]
                    else:
                        formula_sum_dct[_elem] = _p_dct[_elem]

            print(formula_sum_dct)

            merged_idx_lst = []
            # Add Na and K to the sort list
            elem_type_lst = formula_sum_dct.keys()
            for _elem_type in ['C', 'H', 'D', 'O', 'P', 'N', 'Na', 'K', 'S', 'Cl', 'Ca', 'Fe', 'Cu']:
                if _elem_type in elem_type_lst:
                    merged_idx_lst.append(_elem_type)
                else:
                    pass
            # get formula as text
            for _tmp_elem in merged_idx_lst:
                if _tmp_elem in formula_sum_dct.keys():
                    _tmp_elem_num = formula_sum_dct[_tmp_elem]
                    _tmp_info = _tmp_elem + str(_tmp_elem_num)
                    formula_sum_str += _tmp_info
                else:
                    pass

            return formula_sum_str

        else:
            formula_sum_str = 'Error: input is NOT a list'
            return formula_sum_str

    def merge_dct(self, _usr_formula, usr_dct):

        formula_sum_str = ''
        formula_sum_dct = {}

        masscalc = Elem2Mass()

        _parsed_formula = masscalc.get_elem(_usr_formula)

        _p_formula_dct = self.formula2dct(_parsed_formula)
        for _p_dct in [_p_formula_dct, usr_dct]:
            for _elem in _p_dct.keys():
                if _elem in formula_sum_dct.keys():
                    formula_sum_dct[_elem] += _p_dct[_elem]
                else:
                    formula_sum_dct[_elem] = _p_dct[_elem]

        merged_idx_lst = []
        # Add Na and K to the sort list
        elem_type_lst = formula_sum_dct.keys()
        for _elem_type in ['C', 'H', 'D', 'O', 'P', 'N', 'Na', 'K', 'S', 'Cl', 'Ca', 'Fe', 'Cu']:
            if _elem_type in elem_type_lst:
                merged_idx_lst.append(_elem_type)
            else:
                pass
        # get formula as text
        for _tmp_elem in merged_idx_lst:
            if _tmp_elem in formula_sum_dct.keys():
                _tmp_elem_num = formula_sum_dct[_tmp_elem]
                _tmp_info = _tmp_elem + str(_tmp_elem_num)
                formula_sum_str += _tmp_info
            else:
                pass

        return formula_sum_str

    def substract(self, formula_main, formula_part):

        formula_reduced_str = ''

        masscalc = Elem2Mass()
        parsed_formula_lst = []

        _parsed_formula_main = masscalc.get_elem(formula_main)
        _parsed_formula_part = masscalc.get_elem(formula_part)

        formula_main_dct = self.formula2dct(_parsed_formula_main)
        formula_part_dct = self.formula2dct(_parsed_formula_part)

        formula_reduced_dct = copy.deepcopy(formula_main_dct)

        for _elem in formula_part_dct.keys():
            if _elem in formula_main_dct.keys():
                formula_reduced_dct[_elem] -= formula_part_dct[_elem]
            else:
                formula_reduced_str = 'Error: The main formula do not contain: ' + _elem
                return formula_reduced_str

        # remove element with 0 count
        for _elem in formula_reduced_dct.keys():
            if formula_reduced_dct[_elem] == 0:
                del formula_reduced_dct[_elem]

        merged_idx_lst = []
        # Add Na and K to the sort list
        elem_type_lst = formula_reduced_dct.keys()
        for _elem_type in ['C', 'H', 'D', 'O', 'P', 'N', 'Na', 'K', 'S', 'Cl', 'Ca', 'Fe', 'Cu']:
            if _elem_type in elem_type_lst:
                merged_idx_lst.append(_elem_type)
            else:
                pass
        # get formula as text
        for _tmp_elem in merged_idx_lst:
            if _tmp_elem in formula_reduced_dct.keys():
                _tmp_elem_num = formula_reduced_dct[_tmp_elem]
                _tmp_info = _tmp_elem + str(_tmp_elem_num)
                formula_reduced_str += _tmp_info
            else:
                pass

        print('formula_reduced_str', formula_reduced_str)

        return formula_reduced_str

    def substract_dct(self, formula_main, formula_part_dct):

        formula_reduced_str = ''

        masscalc = Elem2Mass()
        parsed_formula_lst = []

        _parsed_formula_main = masscalc.get_elem(formula_main)

        formula_main_dct = self.formula2dct(_parsed_formula_main)

        formula_reduced_dct = copy.deepcopy(formula_main_dct)

        for _elem in formula_part_dct.keys():
            if _elem in formula_main_dct.keys():
                formula_reduced_dct[_elem] -= formula_part_dct[_elem]
            else:
                formula_reduced_str = 'Error: The main formula do not contain: ' + _elem
                return formula_reduced_str

        # remove element with 0 count
        for _elem in formula_reduced_dct.keys():
            if formula_reduced_dct[_elem] == 0:
                del formula_reduced_dct[_elem]

        merged_idx_lst = []
        # Add Na and K to the sort list
        elem_type_lst = formula_reduced_dct.keys()
        for _elem_type in ['C', 'H', 'D', 'O', 'P', 'N', 'Na', 'K', 'S', 'Cl', 'Ca', 'Fe', 'Cu']:
            if _elem_type in elem_type_lst:
                merged_idx_lst.append(_elem_type)
            else:
                pass
        # get formula as text
        for _tmp_elem in merged_idx_lst:
            if _tmp_elem in formula_reduced_dct.keys():
                _tmp_elem_num = formula_reduced_dct[_tmp_elem]
                _tmp_info = _tmp_elem + str(_tmp_elem_num)
                formula_reduced_str += _tmp_info
            else:
                pass

        return formula_reduced_str

# # formula_lst = ['C3H9N', 'C16H36O2Ca', 'C5H14NO4PNa', 'C8H18O3K']
# formula_lst = ['C3H9NO', 'Ca', 'Na', 'PK']
# # formula_lst = []
#
# f = FormulaCalc()
# formula = f.merge(formula_lst)
# formula2 = f.substract(formula_lst[0], 'H2O')
#
# print formula
# print formula2
