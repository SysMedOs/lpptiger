# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re
# import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from DBE2LPP import DBEox


class FAtoLPP(object):

    def __init__(self, usr_smiles):

        fa_smiles_retxt = re.compile(r'(OC\()([CNOPS=\\/@+-]+)(\)=O)')
        fa_smiles_checker = fa_smiles_retxt.match(usr_smiles)
        if fa_smiles_checker:
            self.usr_smiles = usr_smiles
            fa_db_retxt = re.compile(r'(C/C=C\\)')

            self.fa_seg_lst = fa_db_retxt.split(self.usr_smiles)
        else:
            self.usr_smiles = ''
            self.fa_seg_lst = []

        # print 'fa_seg_lst', self.fa_seg_lst

    def get_lpp_all(self):

        fa_seg_lst = self.fa_seg_lst
        print fa_seg_lst
        s_temp_str = ''
        s_pre_oap_lst = []
        s_pre_ocp_lst = []
        s_mid_lst = []

        for _s in fa_seg_lst:
            print 0, '_s', _s
            if _s != 'C/C=C\\' and fa_seg_lst.index(_s) == 0:
                s_temp_str += _s
                s_mid_lst.append(s_temp_str)
                print 0, s_mid_lst

            if _s == 'C/C=C\\' and _s != fa_seg_lst[-1]:
                print 1, '_s', _s
                db_obj = DBEox(_s)
                y = db_obj.get_all_ox(exact=False)
                print 1, y
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

            if _s != 'C/C=C\\' and fa_seg_lst.index(_s) == len(fa_seg_lst) - 1:
                print _s
                for _tmp_oap in s_mid_lst:
                    _s_oap = _tmp_oap + _s
                    s_pre_oap_lst.append(_s_oap)

        print s_pre_oap_lst
        print s_pre_ocp_lst
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

        img = Draw.MolsToGridImage(mol_lst, molsPerRow=6, subImgSize=(300, 300))
        # img = Draw.MolToImage(z_mol, size=(400, 400))
        img.save('New_fa2.png')

        print len(mol_lst), 'image saved!'

s = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'

f_obj = FAtoLPP(s)

f_obj.get_lpp_all()

