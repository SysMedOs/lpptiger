# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of TheoLPP.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

import re
# import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import natsort
from DBEtoLPP import DBEox, LastFAseg


class FAtoLPP(object):

    def __init__(self, usr_smiles):

        fa_smiles_retxt = re.compile(r'(OC\()([CNOPS=\\/@+-]+)(\)=O)')
        fa_smiles_checker = fa_smiles_retxt.match(usr_smiles)
        if fa_smiles_checker:
            self.usr_smiles = usr_smiles
            fa_db_retxt = re.compile(r'(C/C=C\\)')

            self.fa_seg_lst = fa_db_retxt.split(self.usr_smiles)
            # self.fa_seg_lst.remove('')
        else:
            self.usr_smiles = ''
            self.fa_seg_lst = []

    def get_lpp_all(self):

        fa_seg_lst = self.fa_seg_lst
        print 'Split C=C:', fa_seg_lst
        s_temp_str = ''
        s_pre_oap_lst = []
        s_pre_ocp_lst = []
        s_mid_lst = []

        mol_lst = []
        smiles_lst = []
        d_lst = []

        # Define FA acid pattern
        hg = 'OC(CC)=O'
        hg_mol = Chem.MolFromSmiles(hg)
        AllChem.Compute2DCoords(hg_mol)

        if len(fa_seg_lst) > 1:

            for _s in fa_seg_lst:
                # define a dct to sum up all modifications
                # mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0, 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                if _s != 'C/C=C\\' and fa_seg_lst.index(_s) == 0:
                    # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0, 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                    _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                    s_temp_str += _s
                    s_mid_lst.append((s_temp_str, _tmp_mod_dct))

                if _s == 'C/C=C\\' and _s != fa_seg_lst[-1]:
                    db_obj = DBEox(_s)
                    y = db_obj.get_all_ox(exact=False)
                    n_s_mid_lst = []
                    # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0, 'epoxy': 0, 'keto': 0}
                    for _ox_dbe in y:

                        if _ox_dbe[2] == 'OCP':
                            if len(s_mid_lst) == 0:
                                # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0,
                                #                 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                _ox_s_temp_str = s_temp_str + _ox_dbe[0]
                                for _key in _ox_dbe[1].keys():
                                    if _key in _tmp_mod_dct.keys():
                                        _tmp_mod_dct[_key] += _ox_dbe[1][_key]
                                    else:
                                        pass
                                s_pre_ocp_lst.append((_ox_s_temp_str, _tmp_mod_dct))

                            elif len(s_mid_lst) > 0:
                                for _tmp_oap in s_mid_lst:
                                    # _ox_s_temp_str = s_temp_str + _ox_dbe[0]
                                    # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0,
                                    #                 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                    _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                    for _key in _tmp_mod_dct.keys():
                                        if _key in _ox_dbe[1].keys():
                                            _tmp_mod_dct[_key] += _ox_dbe[1][_key]
                                        if _key in _tmp_oap[1].keys():
                                            _tmp_mod_dct[_key] += _tmp_oap[1][_key]
                                        else:
                                            pass
                                    # print '_tmp_oap', _tmp_oap
                                    _ox_s_temp_str = _tmp_oap[0] + _ox_dbe[0]
                                    s_pre_ocp_lst.append((_ox_s_temp_str, _tmp_mod_dct))

                        elif _ox_dbe[2] == 'OAP':

                            for _tmp_oap in s_mid_lst:
                                # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0,
                                #                 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                                for _key in _tmp_mod_dct.keys():
                                    if _key in _ox_dbe[1].keys():
                                        _tmp_mod_dct[_key] += _ox_dbe[1][_key]
                                    if _key in _tmp_oap[1].keys():
                                        _tmp_mod_dct[_key] += _tmp_oap[1][_key]
                                    else:
                                        pass
                                _ox_s_temp_str = _tmp_oap[0] + _ox_dbe[0]
                                n_s_mid_lst.append((_ox_s_temp_str, _tmp_mod_dct))
                                # print n_s_mid_lst
                    s_mid_lst = n_s_mid_lst

                if _s != 'C/C=C\\' and fa_seg_lst.index(_s) == len(fa_seg_lst) - 1:
                    # print _s
                    for _tmp_oap in s_mid_lst:
                        _s_oap = _tmp_oap[0] + _s
                        s_pre_oap_lst.append((_s_oap, _tmp_oap[1]))

                        last_seg_obj = LastFAseg(_s)
                        ox_seg_lst = last_seg_obj.get_all_ox()
                        print _s, 'ox_seg_lst', ox_seg_lst
                        for _seg in ox_seg_lst:
                            _s_oap = _tmp_oap[0] + _seg[0]
                            # _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'OOH': 0, 'epoxy': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                            _tmp_mod_dct = {'dbe': 0, 'OH': 0, 'keto': 0, 'CHO': 0, 'COOH': 0}
                            for _key in _tmp_mod_dct.keys():
                                if _key in _seg[1].keys():
                                    _tmp_mod_dct[_key] += _seg[1][_key]
                                if _key in _tmp_oap[1].keys():
                                    _tmp_mod_dct[_key] += _tmp_oap[1][_key]
                                else:
                                    pass
                            s_pre_oap_lst.append((_s_oap, _tmp_mod_dct))

            s_fin_ocp_lst = []
            for _s in s_pre_ocp_lst:
                s_fin_ocp_lst.append((_s[0] + ')=O', _s[1]))

            s_all_lst = s_fin_ocp_lst + s_pre_oap_lst

            # perform final check
            s_error1 = re.compile(r'C//C')
            rm_lst = []
            for _s_fin in s_all_lst:
                s_checker1 = s_error1.search(_s_fin[0])
                if s_checker1:
                    rm_lst.append(_s_fin)
                else:
                    continue
            for _s_rm in rm_lst:
                s_all_lst.remove(_s_rm)
                print 'removed outlier: ', _s_rm

            for _s in s_all_lst:
                try:
                    _c_counter = str(list(_s[0]).count('C'))
                    _oxfa_str = _c_counter + ':' + str(_s[1]['dbe'])
                    _mod_d_lst = []
                    end_type = ''
                    # for _k in ['OH', 'OOH', 'epoxy', 'keto']:
                    for _k in ['OH', 'keto']:
                        if _s[1][_k] > 0:
                            _mod_d_lst.append(str(_s[1][_k]) + 'x' + _k)
                        else:
                            pass

                    for _e in ['CHO', 'COOH']:
                        if _s[1][_e] > 0:
                            end_type += '(' + _e + '@C' + str(list(_s[0]).count('C')) + ')'
                        else:
                            pass
                    if len(_mod_d_lst) > 0:
                        _mod_str = ','.join(_mod_d_lst)
                        _mod_str = '[' + _mod_str + ']'
                        _oxfa_str += _mod_str
                    else:
                        pass

                    if end_type == '':
                        pass
                    else:
                        _oxfa_str += end_type

                    if _oxfa_str in d_lst:
                        pass

                    else:
                        d_lst.append(_oxfa_str)
                        _mol = Chem.MolFromSmiles(_s[0])
                        AllChem.Compute2DCoords(_mol)
                        AllChem.GenerateDepictionMatching2DStructure(_mol, hg_mol)
                        mol_lst.append(_mol)
                        smiles_lst.append(_s[0])

                except SyntaxError:
                    pass

        elif len(fa_seg_lst) == 1:
            d_lst.append(str(list(fa_seg_lst[0]).count('C')) + ':0')
            _mol = Chem.MolFromSmiles(fa_seg_lst[0])
            AllChem.Compute2DCoords(_mol)
            AllChem.GenerateDepictionMatching2DStructure(_mol, hg_mol)
            mol_lst.append(_mol)
            smiles_lst.append(fa_seg_lst[0])

        sum_lst = zip(d_lst, mol_lst, smiles_lst)
        natsort.natsorted(sum_lst, key=lambda t: t[0])

        # n_mol_lst = []
        # n_d_lst = []
        # for _t in sum_lst:
        #     n_d_lst.append(_t[0])
        #     n_mol_lst.append(_t[1])

        # img = Draw.MolsToGridImage(n_mol_lst, molsPerRow=6, legends=n_d_lst, subImgSize=(300, 300))
        # img.save('New_fa5.png')

        print len(mol_lst), 'image saved!'

        return sum_lst

# s = r'OC(CCCCCCC/C=C\C/C=C\CCCCC)=O'
# # s = r'OC(CCCCCCCCCCCCCCC)=O'
#
# f_obj = FAtoLPP(s)
#
# f_obj.get_lpp_all()

