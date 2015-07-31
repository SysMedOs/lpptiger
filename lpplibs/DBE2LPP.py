# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

# import re
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
            _dbe_hydro = ('C1C(O1)C', {'dbe': 1, '-O-': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C1C(O1)C', {'dbe': 1, '-O-': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CC1C(O1)', {'dbe': 1, '-O-': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_keto(self, exact=False):

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C(=O)CC', {'dbe': 1, '=O': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        elif exact is True and self.dbe is True:
            _dbe_hydro = ('C(=O)CC', {'dbe': 1, '=O': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
            _dbe_hydro = ('CCC(=O)', {'dbe': 1, '=O': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)

        return ox_dbe_lst

    def get_clevage(self):

        return [('C=O', {}, 'OCP'), ('C(O)=O', {}, 'OCP'), ('CC=O', {}, 'OCP'), ('CC(O)=O', {}, 'OCP')]

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
