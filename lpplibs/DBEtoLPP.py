# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem import AllChem, Draw


class DBEox(object):
    """
    Generate all possible lipid peroxidation products from a C=C double bond

    """

    def __init__(self, usr_dbe):
        """
        check if input is a C=C
        :param usr_dbe: str, input should be 'C/C=C\\'
        """
        if usr_dbe == 'C/C=C\\':
            self.dbe = True
        else:
            self.dbe = False

    def get_shift(self):
        """
        shift the 'C/C=C\\' to '/C=C/C'
        'OAP' means the structure will continue extending.
        :return:
        """

        l = [('/C=C/C', {'dbe': 1}, 'OAP')]
        return l

    def get_hydroxy(self, exact=False):

        """
        Generate all hydoxy addition to the structure
        :param exact: False/ True. to use exact position of -OH groups or NOT
                    The use of Exact position will dramatically increase the total number of LPP
        :return: the result list
        """

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C(O)/C=C/', {'dbe': 1, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            # +18 water addition
            # _dbe_hydro = ('C(O)CC', {'dbe': 0, 'OH': 1}, 'OAP')
            # ox_dbe_lst.append(_dbe_hydro)
            # _dbe_hydro = ('C(O)C(O)C', {'dbe': 0, 'OH': 2}, 'OAP')
            # ox_dbe_lst.append(_dbe_hydro)
        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C(O)/C=C/', {'dbe': 1, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('/C=C/C(O)', {'dbe': 1, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('C(O)CC', {'dbe': 0, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CC(O)C', {'dbe': 0, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CCC(O)', {'dbe': 0, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('C(O)C(O)C', {'dbe': 0, 'OH': 2}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CC(O)C(O)', {'dbe': 0, 'OH': 2}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_hydroperoxy(self, exact=False):

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C(OO)/C=C/', {'dbe': 1, 'OOH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C(OO)/C=C/', {'dbe': 1, 'OOH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('/C=C/C(OO)', {'dbe': 1, 'OOH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_epoxy(self, exact=False):

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C1C(O1)C', {'dbe': 0, 'epoxy': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C1C(O1)C', {'dbe': 0, 'epoxy': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CC1C(O1)', {'dbe': 0, 'epoxy': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_keto(self, exact=False):

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C(=O)/C=C/', {'dbe': 1, 'keto': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('C(=O)CC', {'dbe': 0, 'keto': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C(=O)/C=C/', {'dbe': 1, 'keto': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('C(=O)CC', {'dbe': 0, 'keto': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CCC(=O)', {'dbe': 0, 'keto': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_clevage(self):

        return [('C=O', {'CHO': 1}, 'OCP'),
                ('C(O)=O', {'COOH': 1}, 'OCP'),
                ('CC=O', {'CHO': 1}, 'OCP'),
                ('CC(O)=O', {'COOH': 1}, 'OCP')]

    def get_all_ox(self, exact=False):

        if not exact:
            _exact = False
        else:
            _exact = True

        ox_dbe_lst = []
        if self.dbe is True:
            ox_dbe_lst.extend(self.get_hydroxy(exact=_exact))
            ox_dbe_lst.extend(self.get_hydroperoxy(exact=_exact))
            ox_dbe_lst.extend(self.get_epoxy(exact=_exact))
            ox_dbe_lst.extend(self.get_keto(exact=_exact))
            ox_dbe_lst.extend(self.get_clevage())
            ox_dbe_lst.extend(self.get_shift())

        return ox_dbe_lst

class LastFAseg(object):

    def __init__(self, usr_smiles):

        # self.usr_smiles = usr_smiles

        fa_std_code = re.compile(r'(C)([C]+\)=O)')
        fa_smiles_checker = fa_std_code.match(usr_smiles)

        if fa_smiles_checker:
            self.usr_smiles = usr_smiles
            fa_db_retxt = re.compile(r'(C)')

            self.fa_seg_lst = fa_db_retxt.split(self.usr_smiles)
            # self.fa_seg_lst[0] = ''
            # e.g. ['', 'C', '', 'C', '', 'C', '', 'C', '', 'C', ')=O']
            # self.fa_seg_lst.remove('')
        else:
            self.usr_smiles = ''
            self.fa_seg_lst = []

        # print 0, self.fa_seg_lst

    def get_keto(self):

        self.fa_seg_lst[1] = 'C(=O)'
        _ox_seg_str = ''.join(self.fa_seg_lst)
        _ox_seg_lst = [(_ox_seg_str, {'dbe': 0, 'keto': 1}, 'OAP')]

        return _ox_seg_lst

    def get_hydroxy(self):

        self.fa_seg_lst[1] = 'C(O)'
        _ox_seg_str = ''.join(self.fa_seg_lst)
        _ox_seg_lst = [(_ox_seg_str, {'dbe': 0, 'OH': 1}, 'OAP')]

        return _ox_seg_lst

    def get_hydroperoxy(self):

        self.fa_seg_lst[1] = 'C(OO)'
        _ox_seg_str = ''.join(self.fa_seg_lst)
        _ox_seg_lst = [(_ox_seg_str, {'dbe': 0, 'OOH': 1}, 'OAP')]

        return _ox_seg_lst

    def get_all_ox(self):

        _ox_seg_lst = []
        _ox_seg_lst.extend(self.get_hydroxy())
        _ox_seg_lst.extend(self.get_hydroperoxy())
        _ox_seg_lst.extend(self.get_keto())

        # print self.fa_seg_lst
        # print _ox_seg_lst

        return _ox_seg_lst
