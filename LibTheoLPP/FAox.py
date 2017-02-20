# -*- coding: utf-8 -*-
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of TheoLPP.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de

from __future__ import print_function

import re

from temp.PLparser import FAabbr


class FAox(FAabbr):
    def __init__(self, fa_str):
        FAabbr.__init__(self, fa_str)

    def gen_ox_smiles(self, fa_smiles, mode=None):

        fa_smiles_retxt = re.compile(r'(OC\()([CNOPS=\\/@+-]+)(\)=O)')
        fa_smiles_checker = fa_smiles_retxt.match(fa_smiles)
        if fa_smiles_checker:
            fa_slimes_elem = fa_smiles_checker.groups()
            fa_smiles_main = fa_slimes_elem[1]
            print('fa_smiles_main', fa_smiles_main)
            fa_db_retxt = re.compile(r'/C=C\\')

            fa_seg_lst = fa_db_retxt.split(fa_smiles_main)
            tmp_smiles_lst = []
            if len(fa_seg_lst) > 1:
                print('fa_seg_lst', fa_seg_lst)
                tmp_seg_txt = ''
                tmp_seg_ox1_txt = ''
                tmp_seg_ox2_txt = ''
                hydro_counter_0 = 0
                hydro_counter_1 = 0
                hydro_counter_2 = 0
                for tmp_seg in fa_seg_lst:
                    if len(tmp_smiles_lst) == 0:
                        tmp_seg_txt += tmp_seg
                        tmp_smiles_lst.append((tmp_seg_txt, hydro_counter_0))
                        tmp_seg_ox1_txt += tmp_seg
                        tmp_seg_ox2_txt += tmp_seg
                    else:
                        if tmp_seg != fa_seg_lst[-1]:
                            tmp_seg_txt += '/C=C/' + tmp_seg
                            hydro_counter_0 += 0
                            tmp_smiles_lst.append((tmp_seg_txt, hydro_counter_0))
                            if hydro_counter_1 <= 3:
                                tmp_seg_ox1_txt += 'C(O)C' + tmp_seg
                                hydro_counter_1 += 1
                                tmp_smiles_lst.append((tmp_seg_ox1_txt, hydro_counter_1))
                            else:
                                pass
                            if hydro_counter_2 < 2:
                                hydro_counter_2 += 2
                                tmp_seg_ox2_txt += 'C(O)C(O)' + tmp_seg
                                tmp_smiles_lst.append((tmp_seg_ox2_txt, hydro_counter_2))
                            # if hydro_counter_2 == 2:
                            #     hydro_counter_2 += 1
                            #     tmp_seg_ox2_txt += 'C(O)C' + tmp_seg
                            #     tmp_smiles_lst.append((tmp_seg_ox2_txt, hydro_counter_2))
            else:
                print('No DB of this FA!')

            print('tmp_smiles_lst', tmp_smiles_lst)

            smiles_lst = []
            subtxt_lst = []
            for tmp_cmp in tmp_smiles_lst:
                tmp_s = tmp_cmp[0]
                tmp_hydro_counter = tmp_cmp[1]
                if tmp_hydro_counter == 0:
                    tmp_hydro_txt = ''
                else:
                    tmp_hydro_txt = '(' + str(tmp_hydro_counter) + 'x OH)'

                # s_ for SMILES, d_ for Descriptions
                s = 'OC(' + tmp_s + ')=O'

                s_alde1 = 'OC(' + tmp_s + '=O)=O'
                s_alde2 = 'OC(' + tmp_s + 'C=O)=O'
                s_alde3 = 'OC(' + tmp_s + 'CC=O)=O'
                s_acid1 = 'OC(' + tmp_s + '(O)=O)=O'
                s_acid2 = 'OC(' + tmp_s + 'C(O)=O)=O'
                s_acid3 = 'OC(' + tmp_s + 'CC(O)=O)=O'

                s_alde1_count_c = str(list(s_alde1).count('C'))
                s_alde2_count_c = str(list(s_alde2).count('C'))
                s_alde3_count_c = str(list(s_alde3).count('C'))
                s_acid1_count_c = str(list(s_acid1).count('C'))
                s_acid2_count_c = str(list(s_acid2).count('C'))
                s_acid3_count_c = str(list(s_acid3).count('C'))

                mode = 'PL'
                _prefix = ''
                if mode == 'FA':
                    _prefix = 'FA: '
                if mode == 'PL':
                    _prefix = ''
                else:
                    _prefix = ''

                d_alde1_txt = (_prefix + s_alde1_count_c + ':' +
                               str(list(s_alde1).count('=') - 2) + tmp_hydro_txt +
                               '[CHO @C' + s_alde1_count_c + ']')
                d_alde2_txt = (_prefix + s_alde2_count_c + ':' +
                               str(list(s_alde2).count('=') - 2) + tmp_hydro_txt +
                               '[CHO @C' + s_alde2_count_c + ']')
                d_alde3_txt = (_prefix + s_alde3_count_c + ':' +
                               str(list(s_alde3).count('=') - 2) + tmp_hydro_txt +
                               '[CHO @C' + s_alde3_count_c + ']')
                d_acid1_txt = (_prefix + s_acid1_count_c + ':' +
                               str(list(s_acid1).count('=') - 2) + tmp_hydro_txt +
                               '[COOH@C' + s_acid1_count_c + ']')
                d_acid2_txt = (_prefix + s_acid2_count_c + ':' +
                               str(list(s_acid2).count('=') - 2) + tmp_hydro_txt +
                               '[COOH@C' + s_acid2_count_c + ']')
                d_acid3_txt = (_prefix + s_acid3_count_c + ':' +
                               str(list(s_acid3).count('=') - 2) + tmp_hydro_txt +
                               '[COOH@C' + s_acid3_count_c + ']')

                # smiles_lst.append(s)
                smiles_lst.append(s_alde1)
                smiles_lst.append(s_alde2)
                smiles_lst.append(s_alde3)
                smiles_lst.append(s_acid1)
                smiles_lst.append(s_acid2)
                smiles_lst.append(s_acid3)
                subtxt_lst.append(d_alde1_txt)
                subtxt_lst.append(d_alde2_txt)
                subtxt_lst.append(d_alde3_txt)
                subtxt_lst.append(d_acid1_txt)
                subtxt_lst.append(d_acid2_txt)
                subtxt_lst.append(d_acid3_txt)

            # add C(O)/C=C/ species to the list
            _fa_c_db_retxt = re.compile(r'(CC[(]O[)]CC)')
            _fa_d_retxt = re.compile(r'(\d{1,2})')
            _smiles_lst = smiles_lst[:]
            # idx_shift = 0
            for tmp_s_fa in _smiles_lst:
                _fa_c_seg_lst = []
                _fa_c_seg_lst = _fa_c_db_retxt.split(tmp_s_fa)
                # print _fa_c_seg_lst
                _fa_cmod_seg_lst = _fa_c_seg_lst[:]
                for _seg in _fa_cmod_seg_lst:
                    other_ox_keto = 0
                    if _seg == 'CC(O)CC':
                        tmp_idx = _fa_c_seg_lst.index(_seg)
                        other_ox_keto += 1
                        for other_ox in [('C(O)/C=C/C', 1, 'OH'), ('C(=O)/C=C/C', 1, 'C=O'), ('C(=O)CCC', 0, 'C=O')]:

                            _fa_cmod_seg_lst[tmp_idx] = other_ox[0]
                            tmp_s = ''.join(_fa_cmod_seg_lst)
                            tmp_main_idx = smiles_lst.index(tmp_s_fa)
                            # print 'tmp_s', tmp_s
                            smiles_lst.insert(tmp_main_idx, tmp_s)
                            tmp_s_count_c = str(list(tmp_s).count('C'))
                            _tmp_d = subtxt_lst[tmp_main_idx]
                            tmp_d_lst = _fa_d_retxt.split(_tmp_d)
                            # tmp_d_lst ['', '11', ':', '0', '(', '1', 'x OH)[CHO @C', '11', ']']

                            tmp_db = int(tmp_d_lst[3])
                            if tmp_db >= 0:
                                tmp_d_lst[3] = str(tmp_db + other_ox[1])

                            if other_ox[2] == 'OH':
                                pass
                            if other_ox[2] == 'C=O':
                                tmp_hydro = int(tmp_d_lst[5])
                                print(tmp_hydro, 'tmp_d_lst', tmp_d_lst)
                                if tmp_d_lst[6] == 'x OH)[CHO @C' and tmp_hydro > 1:
                                    tmp_d_lst[5] = str(tmp_db - 1)
                                    tmp_d_lst[6] = 'x OH, ' + str(other_ox_keto) + 'x C=O)[CHO @C'
                                if tmp_d_lst[6] == 'x OH)[CHO @C' and tmp_hydro == 1:
                                    tmp_d_lst[5] = str(other_ox_keto)
                                    tmp_d_lst[6] = 'x C=O)[CHO @C'
                                if tmp_d_lst[6] == 'x OH)[COOH@C' and tmp_hydro > 1:
                                    tmp_d_lst[5] = str(tmp_db - 1)
                                    tmp_d_lst[6] = 'x OH, ' + str(other_ox_keto) + 'x C=O)[COOH @C'
                                if tmp_d_lst[6] == 'x OH)[COOH@C' and tmp_hydro == 1:
                                    tmp_d_lst[5] = str(other_ox_keto)
                                    tmp_d_lst[6] = 'x C=O)[COOH @C'
                                # else:
                                #     if
                                #     tmp_d_lst[4] = '(' + str(other_ox_keto) + 'x C=O)'
                            print('tmp_d_lst_ED', tmp_d_lst)
                            tmp_d = ''.join(tmp_d_lst)
                            subtxt_lst.insert(tmp_main_idx, tmp_d)
                            # idx_shift += 1

                    else:
                        pass

            # add OH for the full FA
            _fa_db_retxt = re.compile(r'(C/C=C\\)')
            _fa_seg_lst = _fa_db_retxt.split(fa_smiles_main)

            _fa_db2_retxt = re.compile(r'(C=C)')
            _fa_keto_retxt = re.compile(r'(C\(=O\))')
            _fa_hydro_retxt = re.compile(r'(C\(O\))')


            tmp_smiles_lst = []
            tmp_seg_lst = _fa_seg_lst[:]
            # tmp2_seg_lst = _fa_seg_lst[:]
            oxadd_seg_lst = [[_fa_seg_lst[:]]]
            full_count_c = str(list(fa_smiles).count('C'))
            # if len(full_count_c) == 1:
            #     full_count_c = '_' + full_count_c
            tmp_hydro_counter = 0
            tmp_keto_counter = 0
            tmp_epoxy_counter = 0
            tmp_peroxy_counter = 0

            if 'C/C=C\\' in tmp_seg_lst:
                if tmp_hydro_counter <= 3:
                    for tmp_seg in tmp_seg_lst:

                        if tmp_seg == 'C/C=C\\':
                            tmp_idx = tmp_seg_lst.index(tmp_seg)
                            tmp_seg_lst[tmp_idx] = 'DB processed'
                            print(tmp_idx)
                            current_oxadd_lst = oxadd_seg_lst[:]
                            tmp_hydro_counter += 1
                            for tmp_oxadd_lst in current_oxadd_lst[-1]:
                                print(tmp_idx, 'tmp_oxadd_lst', tmp_oxadd_lst)
                                print('get DB')
                                new_oxadd_lst = []

                                # +16 OH and +18 water addition
                                for oxadd in ['C(O)/C=C/', 'C(O)CC']:
                                    _tmp_oxadd_lst = tmp_oxadd_lst[:]
                                    _tmp_oxadd_lst[tmp_idx] = oxadd

                                    tmp_hydro_txt = '(' + str(tmp_hydro_counter) + 'x OH)'

                                    _tmp_main_fa_txt = ''.join(_tmp_oxadd_lst)
                                    print(_tmp_main_fa_txt)
                                    _tmp_seg_lst = _fa_db2_retxt.split(_tmp_main_fa_txt)
                                    _db_count = _tmp_seg_lst.count('C=C')
                                    _tmp_seg_lst = _fa_hydro_retxt.split(_tmp_main_fa_txt)
                                    _hydro_count = _tmp_seg_lst.count('C(O)')
                                    print('hydro1', _tmp_seg_lst)
                                    _tmp_seg_lst = _fa_keto_retxt.split(_tmp_main_fa_txt)
                                    _keto_count = _tmp_seg_lst.count('C(=O)')
                                    print('keto1', _tmp_seg_lst)
                                    print('Count1', _db_count, _hydro_count, _keto_count)
                                    _tmp_hydro_txt = ''
                                    if _hydro_count >0 and _keto_count > 0:
                                        _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH, ' +
                                                             str(_keto_count) + 'x C=O)')
                                    if _hydro_count >0 and _keto_count == 0:
                                        _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH)')
                                    if _hydro_count ==0 and _keto_count > 0:
                                        _tmp_hydro_txt = ('(' + str(_keto_count) + 'x C=O)')

                                    tmp_s_fa_txt = 'OC(' + _tmp_main_fa_txt + ')=O'
                                    tmp_d_fa_txt = (_prefix + full_count_c + ':' +
                                                    str(_db_count) + _tmp_hydro_txt)
                                    smiles_lst.append(tmp_s_fa_txt)
                                    subtxt_lst.append(tmp_d_fa_txt)
                                    print(tmp_idx, tmp_hydro_counter, 'oxadd_seg_lst', oxadd_seg_lst)
                                    if _tmp_oxadd_lst in new_oxadd_lst:
                                        print('pass')
                                    else:
                                        new_oxadd_lst.append(_tmp_oxadd_lst)

                                if tmp_hydro_counter >= 1:
                                    tmp_hydro_counter_keto = tmp_hydro_counter - 1
                                    tmp_keto_counter += 1

                                    if tmp_hydro_counter_keto > 0 and tmp_keto_counter > 0:
                                        tmp_hydro_txt = ('(' + str(tmp_hydro_counter) + 'x OH, ' +
                                                         str(tmp_keto_counter) + 'x C=O)')
                                    if tmp_hydro_counter_keto == 0 and tmp_keto_counter > 0:
                                        tmp_hydro_txt = ('(' + str(tmp_keto_counter) + 'x C=O)')

                                    # +14 =O group
                                    # http://lipidlibrary.aocs.org/frying/c-epoxy/index.htm
                                    for oxadd in ['C(=O)/C=C/', 'C(=O)CC']:

                                        _tmp_oxadd_lst = tmp_oxadd_lst[:]
                                        _tmp_oxadd_lst[tmp_idx] = oxadd

                                        _tmp_main_fa_txt = ''.join(_tmp_oxadd_lst)
                                        print(_tmp_main_fa_txt)

                                        _tmp_seg_lst = _fa_db2_retxt.split(_tmp_main_fa_txt)
                                        _db_count = _tmp_seg_lst.count('C=C')
                                        _tmp_seg_lst = _fa_hydro_retxt.split(_tmp_main_fa_txt)
                                        _hydro_count = _tmp_seg_lst.count('C(O)')
                                        print('hydro', _tmp_seg_lst)
                                        _tmp_seg_lst = _fa_keto_retxt.split(_tmp_main_fa_txt)
                                        _keto_count = _tmp_seg_lst.count('C(=O)')
                                        print('keto', _tmp_seg_lst)
                                        print('Count', _db_count, _hydro_count, _keto_count)
                                        _tmp_hydro_txt = ''
                                        if _hydro_count > 0 and _keto_count > 0:
                                            _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH, ' +
                                                             str(_keto_count) + 'x C=O)')
                                        if _hydro_count > 0 and _keto_count == 0:
                                            _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH)')
                                        if _hydro_count == 0 and _keto_count > 0:
                                            _tmp_hydro_txt = ('(' + str(_keto_count) + 'x C=O)')

                                        print('_tmp_hydro_txt', _tmp_hydro_txt)

                                        tmp_s_fa_txt = 'OC(' + _tmp_main_fa_txt + ')=O'
                                        tmp_d_fa_txt = (_prefix + full_count_c + ':' +
                                                        str(_db_count) + _tmp_hydro_txt)
                                        smiles_lst.append(tmp_s_fa_txt)
                                        subtxt_lst.append(tmp_d_fa_txt)
                                        print(tmp_idx, tmp_hydro_counter, 'oxadd_seg_lst', oxadd_seg_lst)
                                        if _tmp_oxadd_lst in new_oxadd_lst:
                                            print('pass')
                                        else:
                                            new_oxadd_lst.append(_tmp_oxadd_lst)
                                    oxadd_seg_lst.append(new_oxadd_lst)

                                # add OH after last C=C double bound
                                if 'C/C=C\\' in tmp_seg_lst[(tmp_idx+1):]:
                                    pass
                                else:
                                    if tmp_hydro_counter < 3:
                                        tmp_hydro_counter2 = tmp_hydro_counter
                                        tmp_hydro_counter2 += 1
                                        for oxadd2 in new_oxadd_lst:
                                            tmp_rest = oxadd2[(tmp_idx+1)]
                                            # +16 OH group
                                            tmp_rest_lst = list(tmp_rest)
                                            tmp_rest_lst.insert(1, '(O)')
                                            tmp_rest2 = ''.join(tmp_rest_lst)
                                            oxadd2[(tmp_idx+1)] = tmp_rest2
                                            tmp2_hydro_txt = '(' + str(tmp_hydro_counter2) + 'x OH)'
                                            _tmp_full_fa2_txt = ''.join(oxadd2)
                                            _tmp_seg_lst = _fa_db2_retxt.split(_tmp_full_fa2_txt)
                                            _db_count = _tmp_seg_lst.count('C=C')
                                            _tmp_seg_lst = _fa_hydro_retxt.split(_tmp_full_fa2_txt)
                                            _hydro_count = _tmp_seg_lst.count('C(O)')
                                            _tmp_seg_lst = _fa_keto_retxt.split(_tmp_full_fa2_txt)
                                            _keto_count = _tmp_seg_lst.count('C(=O)')
                                            _tmp_hydro_txt = ''
                                            if _hydro_count >0 and _keto_count > 0:
                                                _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH, ' +
                                                                 str(_keto_count) + 'x C=O)')
                                            if _hydro_count >0 and _keto_count == 0:
                                                _tmp_hydro_txt = ('(' + str(_hydro_count) + 'x OH)')
                                            if _hydro_count ==0 and _keto_count > 0:
                                                _tmp_hydro_txt = ('(' + str(_keto_count) + 'x C=O)')

                                            tmp_s_fa2_txt = 'OC(' + _tmp_full_fa2_txt + ')=O'
                                            tmp_d_fa2_txt = (_prefix + full_count_c + ':' +
                                                             str(_db_count) + _tmp_hydro_txt)
                                            smiles_lst.append(tmp_s_fa2_txt)
                                            subtxt_lst.append(tmp_d_fa2_txt)
                                            # # +14 =O group
                                            # # http://lipidlibrary.aocs.org/frying/c-epoxy/index.htm
                                            # tmp_rest_keto_lst = list(tmp_rest)
                                            # tmp_rest_keto_lst.insert(1, '(=O)')
                                            # tmp_rest2_keto = ''.join(tmp_rest_keto_lst)
                                            # oxadd2[(tmp_idx+1)] = tmp_rest2_keto
                                            # tmp2_hydro_txt = '(' + str(tmp_hydro_counter2) + 'x OH)'
                                            # _tmp_full_fa2_txt = ''.join(oxadd2)
                                            # tmp_s_fa2_txt = 'OC(' + _tmp_full_fa2_txt + ')=O'
                                            # tmp_d_fa2_txt = ('FA: ' + full_count_c + ':' +
                                            #                  str(list(tmp_s_fa2_txt).count('=') - 1) + tmp2_hydro_txt)
                                            # smiles_lst.append(tmp_s_fa2_txt)
                                            # subtxt_lst.append(tmp_d_fa2_txt)

                                        # tmp_hydro_counter = tmp_hydro_counter2

            # remove OH before C=O or COOH
            _fa_c_del_retxt = re.compile(r'(C?[(]O[)]C=O)|(C?[(]O[)]?C[(]O[)]=O)|'
                                         r'(C?[(]O[)]?CC=O)|(C?[(]O[)]?CC[(]O[)]=O)|'
                                         r'(C[(]O[)]C[(]O[)]CC=O)|(C[(]O[)]C[(]O[)]CC[(]O[)]=O)')
            _smiles_lst = smiles_lst[:]
            _subtxt_lst = subtxt_lst[:]
            for tmp_s_del in _smiles_lst:

                tmp_main_idx = smiles_lst.index(tmp_s_del)
                tmp_d_del = subtxt_lst[tmp_main_idx]

                if _fa_c_del_retxt.search(tmp_s_del):
                    smiles_lst.remove(tmp_s_del)
                    subtxt_lst.remove(tmp_d_del)
                else:
                    pass

            return smiles_lst, subtxt_lst

# # #
# # # 'OC(CCCCCCC/C=C\\CCCCCCCC)=O'
# # #
# # x = '20:4(5-Z;8-Z;11-Z;14-Z)'
# # x = '18:3(6-Z;9-Z;12-Z)'
# x = '18:2(9-Z;12-Z)'
# # x = '18:1(9-Z)'
# # x = '18:0'
# FA = FAox(x)
#
# y = FA.FAstr2SMILES(exact_pos=1)
#
# print 'y', y
#
# z_mol = Chem.MolFromSmiles(y[0])
# AllChem.Compute2DCoords(z_mol)
#
# # img = Draw.MolToImage(z_mol)
# # img.show()
#
# (smiles_lst, subtxt_lst) = FA.gen_ox_smiles(y[0])
# print 'smiles_lst', smiles_lst
# print 'subtxt_lst', subtxt_lst
#
# p = 'OC(CC)=O'
# p_mol = Chem.MolFromSmiles(p)
# AllChem.Compute2DCoords(p_mol)
# AllChem.GenerateDepictionMatching2DStructure(z_mol, p_mol)
# z_mol.SetProp("_Name", x)
# s_mol_lst = [z_mol]
# # s_mol_lst = []
# for s in smiles_lst:
#     s_idx = smiles_lst.index(s)
#     tmp_d = subtxt_lst[s_idx]
#     s_mol = Chem.MolFromSmiles(s)
#     AllChem.Compute2DCoords(s_mol)
#     AllChem.GenerateDepictionMatching2DStructure(s_mol, p_mol)
#     s_mol.SetProp("_Name", tmp_d)
#     s_mol_lst.append(s_mol)
#
#
# legends = [x]
# legends += subtxt_lst
# row_num = len(s_mol_lst)//4
# if 1 <= row_num <= 9:
#     pass
# else:
#     row_num = 9
# img2 = Draw.MolsToGridImage(s_mol_lst, molsPerRow=row_num, legends=legends)
# img2.show()
#
# print Chem.MolToMolBlock(s_mol_lst[-1])
#
# w = Chem.SDWriter('fa_test.sdf')
# for m in s_mol_lst:
#     w.write(m)
#
# w.close()
#
# print 'SDF created!'
