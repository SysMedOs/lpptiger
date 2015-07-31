# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
# import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


class DBE(object):

    def __init__(self, usr_dbe):
        if usr_dbe == 'C/C=C\\':
            self.dbe = True
        else:
            self.dbe = False

    def get_shift(self):

        l = [('/C=C/C', {'dbe': 1}, 'OAP')]
        return l

    def get_hydroxy(self, exact=False):

        if not exact:
            exact = False
        else:
            pass

        ox_dbe_lst = []

        if exact is False and self.dbe is True:
            _dbe_hydro = ('C(O)/C=C/', {'dbe': 1, 'OH': 1}, 'OAP')
            ox_dbe_lst.append(_dbe_hydro)
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
            ox_dbe_lst.extend(self.get_clevage())
            ox_dbe_lst.extend(self.get_shift())

        return ox_dbe_lst

s = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
x = r'(C/C=C\\)'

dbe_checker = re.compile(x)

s_lst = dbe_checker.split(s)

print s_lst
mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0, '-O-': 0}
s_pre_oap_lst = []
s_pre_ocp_lst = []
s_mid_lst = []
s_temp_str = ''
for _s in s_lst:
    if _s != 'C/C=C\\' and s_lst.index(_s) == 0:
        s_temp_str += _s
        s_mid_lst.append(s_temp_str)
        print 0, s_mid_lst

    if _s == 'C/C=C\\' and _s != s_lst[-1]:
        print _s
        db_obj = DBE(_s)
        y = db_obj.get_all_ox(exact=False)
        # print 1, y
        n_s_mid_lst = []
        for _ox_dbe in y:
            if _ox_dbe[2] == 'OCP':
                if len(s_mid_lst) == 0:
                    _ox_s_temp_str = s_temp_str + _ox_dbe[0]
                    s_pre_ocp_lst.append(_ox_s_temp_str)

                elif len(s_mid_lst) > 0:
                    for _tmp_oap in s_mid_lst:
                        # print '_tmp_oap', _tmp_oap
                        _ox_s_temp_str = _tmp_oap + _ox_dbe[0]
                        s_pre_ocp_lst.append(_ox_s_temp_str)

            elif _ox_dbe[2] == 'OAP':

                for _tmp_oap in s_mid_lst:
                    _ox_s_temp_str = _tmp_oap + _ox_dbe[0]
                    n_s_mid_lst.append(_ox_s_temp_str)
                    # print n_s_mid_lst
        s_mid_lst = n_s_mid_lst
        print 's_mid_lst', s_mid_lst

    if _s != 'C/C=C\\' and s_lst.index(_s) == len(s_lst) - 1:
        print _s
        for _tmp_oap in s_mid_lst:
            _s_oap = _tmp_oap + _s
            s_pre_oap_lst.append(_s_oap)

s_fin_ocp_lst = []
for _s in s_pre_ocp_lst:
    s_fin_ocp_lst.append(_s + ')=O')
    # s_fin_ocp_lst.append(_s + ')=O')
    # _ox_s_temp_str = _tmp_oap + _s
    # s_pre_oap_lst.append(_ox_s_temp_str)

print 8, s_fin_ocp_lst

print 9, s_pre_oap_lst

s_all_lst = s_fin_ocp_lst + s_pre_oap_lst

hg = 'OC(CCCC)=O'
hg_mol = Chem.MolFromSmiles(hg)
AllChem.Compute2DCoords(hg_mol)

mol_lst = []
for _s in s_all_lst:
    try:
        _mol = Chem.MolFromSmiles(_s)
        AllChem.Compute2DCoords(_mol)
        AllChem.GenerateDepictionMatching2DStructure(_mol, hg_mol)
        mol_lst.append(_mol)
    except:
        pass

img = Draw.MolsToGridImage(mol_lst, molsPerRow=7, subImgSize=(300, 300))
# img = Draw.MolToImage(z_mol, size=(400, 400))
img.save('New_fa.png')

print len(mol_lst), 'image saved!'
