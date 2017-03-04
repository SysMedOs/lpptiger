# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPsmi.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from ExactMassCalc import Elem2Mass

class SMILESparser(object):

    def __init__(self):

        self.charge_dct = {'M': {}, '[M+H]+': {'H': 1}, '[M+Na]+': {'Na': 1}, '[M+K]+': {'K': 1},
                           '[M-H2O+H]+': {'H': -1, 'O': -1},
                           '[M-H2O+Na]+': {'H': -2, 'O': -1, 'Na': 1},
                           '[M-H2O+K]+': {'H': -2, 'O': -1, 'K': 1},
                           '[M-H]-': {'H': -1},
                           '[M+FA-H]-': {'C': 1, 'H': 1, 'O': 2},
                           '[M+HCOO]-': {'C': 1, 'H': 1, 'O': 2},
                           '[M-H2O-H]-': {'H': -3, 'O': -1},
                           '[M-H2O+FA-H]-': {'C': 1, 'H': -1, 'O': 1},
                           '[M-H2O+HCOO]-': {'C': 1, 'H': -1, 'O': 1}}
        self.charge_lst = ['M', '[M+H]+', '[M+Na]+', '[M+K]+', '[M-H2O+H]+', '[M-H2O+Na]+', '[M-H2O+K]+',
                           '[M-H]-', '[M+FA-H]-', '[M+HCOO]-', '[M-H2O-H]-', '[M-H2O+FA-H]-', '[M-H2O+HCOO]-']

    def smiles2elem(self, usr_smiles):

        """
        Generate a dict of element : atom_counts
        For Neutral compound only.
        !!! DO NOT support Na, K charged species !!!
        :return: elem_dct
        """

        smiles_lst = list(usr_smiles)
        elem_dct = {}
        elem_dct['C'] = smiles_lst.count('C')
        elem_dct['O'] = smiles_lst.count('O')
        if smiles_lst.count('P') > 0:
            elem_dct['P'] = smiles_lst.count('P')
        if smiles_lst.count('N') > 0:
            elem_dct['N'] = smiles_lst.count('N')
        elem_dct['dbe'] = smiles_lst.count('=')
        elem_dct['H'] = smiles_lst.count('C') * 2 + 2 + 4 - 2 * smiles_lst.count('=')

        return elem_dct

    def smiles2formula(self, usr_smiles, charge=None):

        formula = ''

        if not charge:
            charge = 'M'

        if charge in self.charge_lst:
            pass
        else:
            formula = 'Charge Error: ' + charge
            return formula

        smiles_lst = list(usr_smiles)
        elem_dct = {'C': smiles_lst.count('C'), 'O': smiles_lst.count('O'), 'dbe': smiles_lst.count('=')}
        if smiles_lst.count('P') > 0:
            elem_dct['P'] = smiles_lst.count('P')
        if smiles_lst.count('N') > 0:
            elem_dct['N'] = smiles_lst.count('N')
            elem_dct['H'] = smiles_lst.count('C') * 2 + 2 + 4 - 2 * smiles_lst.count('=') - 9
            elem_dct['D'] = 9
        elif smiles_lst.count('N') == 0:
            elem_dct['H'] = smiles_lst.count('C') * 2 + 2 - 2 * smiles_lst.count('=')

        if charge in self.charge_lst:
            tmp_chg_dct = self.charge_dct[charge]

            for _atom in tmp_chg_dct.keys():
                if _atom in elem_dct.keys():
                    elem_dct[_atom] += tmp_chg_dct[_atom]
                else:
                    elem_dct[_atom] = tmp_chg_dct[_atom]
        if charge == '':
            pass

        _charged_str = ''
        _charged_idx_lst = ['C', 'H', 'D', 'O']
        # Add Na and K to the sort list
        if 'P' in elem_dct.keys():
            _charged_idx_lst.append('P')
        if 'N' in elem_dct.keys():
            _charged_idx_lst.append('N')
        if 'Na' in elem_dct.keys():
            _charged_idx_lst.append('Na')
        if 'K' in elem_dct.keys():
            _charged_idx_lst.append('K')
        else:
            pass

        # get formula as text
        for _tmp_elem in _charged_idx_lst:
            if _tmp_elem in elem_dct.keys():
                _tmp_elem_num = elem_dct[_tmp_elem]
                _tmp_info = _tmp_elem + str(_tmp_elem_num)
                _charged_str += _tmp_info
            else:
                pass

        return _charged_str
