# -*- coding: utf-8 -*-
# Copyright 2014-2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig                #
# The software is currently  under development and is not ready to be released.         #
# A suitable license will be chosen before the official release of LipidSDFcreator.     #
# For more info please contact: zhixu.ni@uni-leipzig.de                                 #

from ExactMassCalc import Elem2Mass
import pandas as pd
import re
import csv
import math


class FAabbr():
    def __init__(self, fa_str):

        self.usr_FA_str = fa_str

    def FAstr2elem(self):

        usr_FA_lst = self.usr_FA_str.split('(')

        # decode FA main chain and Modifications
        if len(usr_FA_lst) == 1:
            FA_chain_str = usr_FA_lst[0]
            FA_mod_str = ''
        elif len(usr_FA_lst) == 2:
            FA_chain_str = usr_FA_lst[0]
            FA_mod_str = usr_FA_lst[1].strip('()')
        else:
            FA_chain_str = ''
            FA_mod_str = ''

        # process FA Chain

        if FA_chain_str != '':
            # lyso
            if FA_chain_str.lower() == 'lyso':
                FA_typ = 'lyso'
                FA_chain_L_num = 0
                FA_chain_DBE_num = 0
            else:
                try:
                    FA_chain_lst = FA_chain_str.split(':')
                    if len(FA_chain_lst) == 2:
                        FA_typ = 'general'
                        FA_chain_L_num = int(FA_chain_lst[0])
                        FA_chain_DBE_num = int(FA_chain_lst[1])

                    else:
                        FA_typ = 'special'
                        FA_chain_L_num = 0
                        FA_chain_DBE_num = 0

                except:
                    FA_typ = 'special'
                    FA_chain_L_num = 0
                    FA_chain_DBE_num = 0
        else:
            FA_typ = 'missing'
            FA_chain_L_num = 0
            FA_chain_DBE_num = 0

        # print FA_chain_L_num, FA_chain_DBE_num, FA_typ

        # process FA modification

        # define modifications
        # mod_DB_dct = {'H': -2}
        mod_O_dct = {'O': 1}
        mod_OH_dct = {'O': 1}
        mod_epoxy_dct = {'O': 1}
        mod_OOH_dct = {'O': 2}
        mod_oxo_dct = {'H': -2, 'O': 1}
        mod_COOH_dct = {'H': -2, 'O': 2}
        mod_Furyl_dct = {'C': 4, 'H': 2, 'O': 1}
        mod_OCH3_dct = {'C': 1, 'H': 2, 'O': 1}
        mod_sum_dct = {'O': mod_O_dct, 'OH': mod_OH_dct, 'epoxy': mod_epoxy_dct,
                       'OOH': mod_OOH_dct, 'oxo': mod_oxo_dct,
                       'COOH': mod_COOH_dct, 'Furyl': mod_Furyl_dct,
                       'furyl': mod_Furyl_dct, 'OCH3': mod_OCH3_dct
                       }

        # Create dct to store modifications
        FA_mod_elem_dct = {'H': 0, 'O': 0}

        if FA_mod_str == '':
            pass
        else:
            FA_mod_lst = FA_mod_str.split(';')

            for tmp_mod in FA_mod_lst:
                tmp_mod = tmp_mod.strip(' ')

                # checker for mod with place
                mod_p_std = re.compile('\d{1,2}-\D*')
                mod_checker_p = mod_p_std.match(tmp_mod)

                # checker for mod with NO specific position
                mod_g_std = re.compile('(\d{1,2})(\D*)')
                mod_checker_g = mod_g_std.match(tmp_mod)
                if mod_checker_p:
                    # print tmp_mod, 'pass mod idx'
                    tmp_mod_lst = tmp_mod.split('-')
                    tmp_mod_p_num = int(tmp_mod_lst[0])
                    tmp_mod_typ = tmp_mod_lst[1]
                    # put Elem number into dict
                    if tmp_mod_typ in mod_sum_dct.keys():
                        for tmp_elem_k, tmp_elem_num in mod_sum_dct[tmp_mod_typ].items():
                            # print tmp_mod_typ, tmp_elem_k, tmp_elem_num
                            if tmp_elem_k in FA_mod_elem_dct.keys():
                                FA_mod_elem_dct[tmp_elem_k] += tmp_elem_num
                            else:
                                FA_mod_elem_dct[tmp_elem_k] = tmp_elem_num

                elif mod_checker_g:
                    # print tmp_mod, 'pass mod num'
                    tmp_mod_cnt_num = int(mod_checker_g.group(1))
                    tmp_mod_typ = mod_checker_g.group(2)
                    if tmp_mod_typ in mod_sum_dct.keys():
                        for tmp_elem_k, tmp_elem_num in mod_sum_dct[tmp_mod_typ].items():
                            # print tmp_mod_typ, tmp_elem_k, tmp_elem_num
                            if tmp_elem_k in FA_mod_elem_dct.keys():
                                FA_mod_elem_dct[tmp_elem_k] += tmp_elem_num * tmp_mod_cnt_num
                            else:
                                FA_mod_elem_dct[tmp_elem_k] = tmp_elem_num
                else:
                    # print tmp_mod, 'Not MOD'
                    pass

        # print FA_mod_elem_dct

        FA_chain_elem_dct = {}
        FA_chain_elem_checker = 0

        if FA_typ in ['special', 'missing']:
            pass
        elif FA_typ in ['lyso']:
            FA_chain_elem_dct['H'] = 2
            FA_chain_elem_dct['O'] = 1
            FA_chain_elem_checker = 1

        elif FA_typ in ['general']:
            # calc chain without mod
            tmp_FA_C_num = FA_chain_L_num
            tmp_FA_H_num = (tmp_FA_C_num - FA_chain_DBE_num) * 2
            tmp_FA_O_num = 2

            FA_chain_elem_dct['C'] = tmp_FA_C_num
            FA_chain_elem_dct['H'] = tmp_FA_H_num
            FA_chain_elem_dct['O'] = tmp_FA_O_num

            for tmp_elem_k, tmp_elem_num in FA_mod_elem_dct.items():
                # print tmp_mod_typ, tmp_elem_k, tmp_elem_num
                if tmp_elem_k in FA_mod_elem_dct.keys():
                    FA_chain_elem_dct[tmp_elem_k] += tmp_elem_num
                else:
                    FA_chain_elem_dct[tmp_elem_k] = tmp_elem_num

            FA_chain_elem_checker = 1

        else:
            FA_chain_elem_dct = {}

        return FA_chain_elem_dct, FA_chain_elem_checker

    def FAstr2SMILES(self, exact_pos=1):

        if exact_pos == 0:
            FA_ExactPos_checker = 0
        else:
            FA_ExactPos_checker = 1

        usr_FA_lst = self.usr_FA_str.split('(')
        FA_SMILES_str = ''

        # decode FA main chain and Modifications
        if len(usr_FA_lst) == 1:
            FA_chain_str = usr_FA_lst[0]
            FA_mod_str = ''
        elif len(usr_FA_lst) == 2:
            FA_chain_str = usr_FA_lst[0]
            FA_mod_str = usr_FA_lst[1].strip('()')
        else:
            FA_chain_str = ''
            FA_mod_str = ''

        # process FA Chain

        if FA_chain_str != '':
            # lyso
            if FA_chain_str.lower() == 'lyso':
                FA_typ = 'lyso'
                FA_chain_L_num = 0
                FA_chain_DBE_num = 0
                FA_ExactPos_checker = 1
            else:
                try:
                    FA_chain_lst = FA_chain_str.split(':')
                    if len(FA_chain_lst) == 2:
                        FA_typ = 'general'
                        FA_chain_L_num = int(FA_chain_lst[0])
                        FA_chain_DBE_num = int(FA_chain_lst[1])

                    else:
                        FA_typ = 'special'
                        FA_chain_L_num = 0
                        FA_chain_DBE_num = 0

                except:
                    FA_typ = 'special'
                    FA_chain_L_num = 0
                    FA_chain_DBE_num = 0
        else:
            FA_typ = 'missing'
            FA_chain_L_num = 0
            FA_chain_DBE_num = 0

        # process FA modification

        # define modifications
        mod_sum_dct = {'DB': r'C=', 'Z': r'C=', 'E': r'C=', 'OH': r'C(O)', 'epoxy': r'C(O[i])C[i]',
                       'OOH': r'C(OO)', 'oxo': r'C([sub])=O',
                       'COOH': r'C(O)=O'}

        # Create dct to store modifications
        FA_SMILES_lst = ['OC(']
        FA_SMILES_end_lst = [')=O']
        FA_chain_lite_lst = ['C'] * (FA_chain_L_num - 1)

        # Create a list to indicate the modification status of each carbon
        FA_Chain_mod_lst = [0 for i in range(FA_chain_L_num - 1)]

        if FA_mod_str == '':
            if FA_chain_str.lower() == 'lyso':
                FA_SMILES_str = 'lyso'
                FA_ExactPos_checker = 1

            else:
                if FA_chain_L_num == 1 and FA_chain_DBE_num == 0:
                    FA_SMILES_str = 'OC=O'

                elif FA_chain_L_num == 2 and FA_chain_DBE_num == 0:
                    FA_SMILES_str = 'OCC(=O)O'

                elif FA_chain_L_num > 2 and FA_chain_DBE_num == 0:
                    FA_SMILES_lst += FA_chain_lite_lst
                    FA_SMILES_lst += FA_SMILES_end_lst
                    FA_SMILES_str = ('').join(FA_SMILES_lst)
                    FA_ExactPos_checker = 1

                elif FA_chain_L_num > 2 and FA_chain_DBE_num > 0:
                    # put double bounds
                    common_DB_dct = {(16, 1): [9], (18, 1): [9], (18, 2): [9, 12], (18, 3): [9, 12, 15],
                                     (20, 4): [5, 8, 11, 14], (22, 6): [4, 7, 10, 13, 16, 19]}

                    print("Generate DB for NO modi LPP")
                    # for common positions
                    if (FA_chain_L_num, FA_chain_DBE_num) in common_DB_dct.keys():
                        for db_pos in common_DB_dct[(FA_chain_L_num, FA_chain_DBE_num)]:
                            if FA_chain_lite_lst[db_pos - 2] == 'C' and FA_chain_lite_lst[db_pos - 1] == 'C':
                                FA_chain_lite_lst[db_pos - 2] = '/' + mod_sum_dct['DB']
                                FA_chain_lite_lst[db_pos - 1] = 'C' + '\\'
                                FA_Chain_mod_lst[db_pos - 2] = 1
                                FA_Chain_mod_lst[db_pos - 1] = 1
                                FA_ExactPos_checker = 1

                        FA_SMILES_lst += FA_chain_lite_lst
                        FA_SMILES_lst += FA_SMILES_end_lst
                        FA_SMILES_str = ('').join(FA_SMILES_lst)

                    else:
                        start_pos = int(math.floor(FA_chain_L_num/(FA_chain_DBE_num + 1)))
                        if start_pos > 1:
                            pass

                        if start_pos == 1:
                            start_pos = 2

                        tmp_db_pos_lst = list(range(FA_chain_DBE_num))
                        db_pos_lst = []
                        for tmp_pos in tmp_db_pos_lst:
                            pos_idx = tmp_db_pos_lst.index(tmp_pos)
                            if start_pos + 3 * FA_chain_DBE_num <= FA_chain_L_num:
                                db_pos_lst.append(start_pos + 3 * pos_idx)

                            elif start_pos + 2 * FA_chain_DBE_num <= FA_chain_L_num:
                                db_pos_lst.append(start_pos + 2)

                            else:
                                db_pos_lst = []
                                FA_SMILES_str = ''

                        for db_pos in db_pos_lst:
                            if FA_chain_lite_lst[db_pos - 2] == 'C' and FA_chain_lite_lst[db_pos - 1] == 'C':
                                FA_chain_lite_lst[db_pos - 2] = '/' + mod_sum_dct['DB']
                                FA_chain_lite_lst[db_pos - 1] = 'C' + '\\'
                                FA_Chain_mod_lst[db_pos - 2] = 1
                                FA_Chain_mod_lst[db_pos - 1] = 1

                        FA_SMILES_lst += FA_chain_lite_lst
                        FA_SMILES_lst += FA_SMILES_end_lst
                        FA_SMILES_str = ('').join(FA_SMILES_lst)

        else:

            FA_mod_lst = FA_mod_str.split(';')
            print('FA_mod_lst', FA_mod_lst)
            FA_subchain_lst = []

            # define a list to store modification positions
            usr_mod_idx_lst = []

            for tmp_mod in FA_mod_lst:
                tmp_mod = tmp_mod.strip(' ')

                # checker for mod with place
                mod_p_std = re.compile('\d{1,2}-\D*')
                mod_checker_p = mod_p_std.match(tmp_mod)

                # checker for mod with NO specific position
                mod_g_std = re.compile('(\d{1,2})(\D*)')
                mod_checker_g = mod_g_std.match(tmp_mod)
                if mod_checker_p and FA_chain_lite_lst != []:
                    # print tmp_mod, 'pass mod idx'
                    tmp_mod_lst = tmp_mod.split('-')
                    tmp_mod_idx = int(tmp_mod_lst[0])
                    tmp_mod_typ = tmp_mod_lst[1]
                    # put Elem number into dict
                    if tmp_mod_idx <= FA_chain_L_num and tmp_mod_typ in mod_sum_dct.keys():
                        usr_mod_idx_lst.append(tmp_mod_idx)
                        if tmp_mod_typ in ['OH', 'OOH']:
                            if FA_chain_lite_lst[tmp_mod_idx - 2] == 'C':
                                FA_chain_lite_lst[tmp_mod_idx - 2] = mod_sum_dct[tmp_mod_typ]
                                FA_Chain_mod_lst[tmp_mod_idx - 2] = 1

                            else:
                                FA_SMILES_lst = []
                                FA_chain_lite_lst = []
                                FA_subchain_lst = []
                                FA_SMILES_str = ''

                        elif tmp_mod_typ in ['DB', 'Z']:
                            if FA_chain_lite_lst[tmp_mod_idx - 2] == 'C' and FA_chain_lite_lst[tmp_mod_idx - 1] == 'C':
                                FA_chain_lite_lst[tmp_mod_idx - 2] = '/' + mod_sum_dct[tmp_mod_typ]
                                FA_chain_lite_lst[tmp_mod_idx - 1] = 'C' + '\\'
                                FA_Chain_mod_lst[tmp_mod_idx - 2] = 1
                                FA_Chain_mod_lst[tmp_mod_idx - 1] = 1

                        elif tmp_mod_typ in ['E']:
                            if FA_chain_lite_lst[tmp_mod_idx - 2] == 'C' and FA_chain_lite_lst[tmp_mod_idx - 1] == 'C':
                                FA_chain_lite_lst[tmp_mod_idx - 2] = '/' + mod_sum_dct[tmp_mod_typ]
                                FA_chain_lite_lst[tmp_mod_idx - 1] = 'C' + '/'
                                FA_Chain_mod_lst[tmp_mod_idx - 2] = 1
                                FA_Chain_mod_lst[tmp_mod_idx - 1] = 1

                        elif tmp_mod_typ == 'oxo':
                            if tmp_mod_idx == FA_chain_L_num:
                                FA_chain_lite_lst[-1] = 'C([H])=O'
                                FA_Chain_mod_lst[-1] = 1

                            else:
                                FA_chain_lite_lst[tmp_mod_idx - 2] = 'C('
                                FA_subchain_lst.append(tmp_mod_idx - 1)
                                FA_Chain_mod_lst[tmp_mod_idx - 2] = 1

                        print(FA_chain_lite_lst)
                        print('FA_Chain_mod_lst', FA_Chain_mod_lst)

                    else:
                        pass

                elif mod_checker_g:
                    if exact_pos == 1:
                        FA_ExactPos_checker = 1 * FA_ExactPos_checker
                        FA_SMILES_lst = []
                        FA_chain_lite_lst = []
                        FA_subchain_lst = []
                        FA_SMILES_str = ''

                    elif exact_pos == 0:
                        FA_ExactPos_checker = 0 * FA_ExactPos_checker
                        print('Try to get structure --->', self.usr_FA_str)

                        # put double bounds
                        common_DB_dct = {(18, 1): [9], (18, 2): [9, 12], (18, 3): [9, 12, 15],
                                         (20, 4): [5, 8, 11, 14], (22, 6): [4, 7, 10, 13, 16, 19]}
                        if FA_chain_DBE_num == 0:
                            pass

                        else:
                            # for common positions
                            if (FA_chain_L_num, FA_chain_DBE_num) in common_DB_dct.keys():
                                for db_pos in common_DB_dct[(FA_chain_L_num, FA_chain_DBE_num)]:
                                    if FA_chain_lite_lst[db_pos - 1] == 'C' and FA_chain_lite_lst[db_pos] == 'C':
                                        FA_chain_lite_lst[db_pos - 1] = '/' + mod_sum_dct['DB']
                                        FA_chain_lite_lst[db_pos] = 'C' + '\\'
                                        FA_Chain_mod_lst[db_pos - 1] = 1
                                        FA_Chain_mod_lst[db_pos] = 1
                            else:
                                start_pos = int(math.floor(FA_chain_L_num/(FA_chain_DBE_num + 1)))
                                if start_pos > 1:
                                    pass

                                if start_pos == 1:
                                    start_pos = 2

                                tmp_db_pos_lst = list(range(FA_chain_DBE_num))
                                db_pos_lst = []
                                for tmp_pos in tmp_db_pos_lst:
                                    pos_idx = tmp_db_pos_lst.index(tmp_pos)
                                    if start_pos + 3 * FA_chain_DBE_num <= FA_chain_L_num:
                                        db_pos_lst.append(start_pos + 3 * pos_idx)

                                    elif start_pos + 2 * FA_chain_DBE_num <= FA_chain_L_num:
                                        db_pos_lst.append(start_pos + 2)

                                    else:
                                        db_pos_lst = []
                                        FA_chain_lite_lst = []
                                        FA_Chain_mod_lst = []
                                        FA_SMILES_lst = []
                                        FA_subchain_lst = []
                                        FA_SMILES_str = ''

                                for db_pos in db_pos_lst:
                                    if FA_chain_lite_lst[db_pos - 2] == 'C' and FA_chain_lite_lst[db_pos - 1] == 'C':
                                        FA_chain_lite_lst[db_pos - 2] = '/' + mod_sum_dct['DB']
                                        FA_chain_lite_lst[db_pos - 1] = 'C' + '\\'
                                        FA_Chain_mod_lst[db_pos - 2] = 1
                                        FA_Chain_mod_lst[db_pos - 1] = 1

                        tmp_mod_num = int(mod_checker_g.group(1))
                        tmp_mod_typ = mod_checker_g.group(2)
                        # if tmp_mod_typ == 'oxo':
                        #     tmp_mod_idx = FA_chain_L_num
                        #     tmp_mod_idx_lst = [tmp_mod_idx]
                        #     new_mod_idx_lst = tmp_mod_idx_lst
                        #
                        # else:
                        #     tmp_mod_idx_lst = range(1, int(int(tmp_mod_num) + 1))
                        #     new_mod_idx_lst = tmp_mod_idx_lst
                        #     for tmp_idx in tmp_mod_idx_lst:
                        #         if tmp_idx in usr_mod_idx_lst:
                        #             if tmp_idx <= FA_chain_L_num:
                        #                 new_mod_idx = tmp_idx
                        #                 while new_mod_idx <= FA_chain_L_num:
                        #                     if new_mod_idx in usr_mod_idx_lst:
                        #                         new_mod_idx += 1
                        #                     else:
                        #                         new_mod_idx_lst = [new_mod_idx if x == tmp_idx
                        #                                            else x for x in tmp_mod_idx_lst]
                        #                         break
                        #             else:
                        #                 print('too many mods')
                        #                 new_mod_idx_lst = []
                        #         else:
                        #             pass
                        # print(tmp_mod_idx_lst, new_mod_idx_lst, usr_mod_idx_lst)

                        # put Elem number into dict
                        # for tmp_mod_idx in new_mod_idx_lst:
                        while tmp_mod_num > 0:
                            if FA_Chain_mod_lst != [] and FA_chain_lite_lst != [] and tmp_mod_typ in mod_sum_dct.keys():
                                # usr_mod_idx_lst.append(tmp_mod_idx)
                                # if tmp_mod_typ in ['OH', 'OOH']:
                                #     if FA_chain_lite_lst[tmp_mod_idx - 1] == 'C':
                                #         FA_chain_lite_lst[tmp_mod_idx - 1] = mod_sum_dct[tmp_mod_typ]
                                if tmp_mod_typ in ['OH', 'OOH']:
                                    if FA_chain_L_num > 4:
                                        tmp_FA_Chain_mod_lst = FA_Chain_mod_lst[3:]
                                        try:
                                            tmp_ins_idx = tmp_FA_Chain_mod_lst.index(0)
                                        except ValueError:
                                            tmp_ins_idx = 0

                                        print(FA_chain_lite_lst, tmp_FA_Chain_mod_lst, tmp_ins_idx)

                                        if FA_chain_lite_lst[tmp_ins_idx + 3] == 'C':
                                            FA_chain_lite_lst[tmp_ins_idx + 3] = mod_sum_dct[tmp_mod_typ]
                                            FA_Chain_mod_lst[tmp_ins_idx + 3] = 1

                                        else:
                                            FA_SMILES_lst = []
                                            FA_subchain_lst = []
                                            FA_SMILES_str = ''

                                    else:
                                        tmp_ins_idx = FA_Chain_mod_lst.index(0)
                                        if FA_chain_lite_lst[tmp_ins_idx] == 'C':
                                            FA_chain_lite_lst[tmp_ins_idx] = mod_sum_dct[tmp_mod_typ]
                                            FA_Chain_mod_lst[tmp_ins_idx] = 1

                                        else:
                                            FA_SMILES_lst = []
                                            FA_subchain_lst = []
                                            FA_SMILES_str = ''

                                    tmp_mod_num += -1

                                # elif tmp_mod_typ in ['DB', 'Z']:
                                #     if FA_chain_DBE_num > 0:
                                #         if FA_chain_lite_lst[tmp_mod_idx - 1] == 'C' and FA_chain_lite_lst[tmp_mod_idx - 1] == 'C':
                                #             FA_chain_lite_lst[tmp_mod_idx - 1] = '/' + mod_sum_dct[tmp_mod_typ]
                                #             FA_chain_lite_lst[tmp_mod_idx - 1] = 'C' + '\\'
                                #     else:
                                #         pass
                                #
                                # elif tmp_mod_typ in ['E']:
                                #     if FA_chain_lite_lst[tmp_mod_idx - 1] == 'C' and FA_chain_lite_lst[tmp_mod_idx - 1] == 'C':
                                #         FA_chain_lite_lst[tmp_mod_idx - 1] = '/' + mod_sum_dct[tmp_mod_typ]
                                #         FA_chain_lite_lst[tmp_mod_idx - 1] = 'C' + '/'

                                elif tmp_mod_typ == 'oxo':
                                    if FA_Chain_mod_lst != [] and FA_Chain_mod_lst[-1] == 0:
                                        FA_chain_lite_lst[-1] = 'C([H])=O'
                                        FA_Chain_mod_lst[-1] = 1

                                    else:
                                        if FA_Chain_mod_lst != [] and FA_chain_L_num > 4:
                                            tmp_FA_Chain_mod_lst = FA_Chain_mod_lst[3:]
                                            tmp_ins_idx = tmp_FA_Chain_mod_lst.index(0)
                                            if FA_chain_lite_lst[tmp_ins_idx + 3] == 'C':
                                                FA_chain_lite_lst[tmp_ins_idx + 3] = 'C('
                                                FA_subchain_lst.append(tmp_ins_idx + 3 - 1)
                                                FA_Chain_mod_lst[tmp_ins_idx + 3] = 1

                                        else:
                                            if FA_Chain_mod_lst != []:
                                                tmp_ins_idx = FA_Chain_mod_lst.index(0)
                                                if FA_chain_lite_lst[tmp_ins_idx] == 'C':
                                                    FA_chain_lite_lst[tmp_ins_idx] = 'C('
                                                    FA_subchain_lst.append(tmp_ins_idx - 1)
                                                    FA_Chain_mod_lst[tmp_ins_idx] = 1

                                    tmp_mod_num += -1

                                elif tmp_mod_typ == 'COOH':
                                    if FA_Chain_mod_lst[-1] == 0:
                                        FA_chain_lite_lst[-1] = 'C([OH])=O'
                                        FA_Chain_mod_lst[-1] = 1
                                    else:
                                        FA_SMILES_lst = []
                                        FA_subchain_lst = []
                                        FA_SMILES_str = ''

                                    tmp_mod_num += -1

                                # print('FA_chain_lite_lst', FA_chain_lite_lst)
                                else:
                                    FA_SMILES_lst = []
                                    FA_subchain_lst = []
                                    FA_SMILES_str = ''
                                    tmp_mod_num += -1
                            else:
                                FA_SMILES_lst = []
                                FA_subchain_lst = []
                                FA_SMILES_str = ''
                                tmp_mod_num += -1
                else:
                    # print tmp_mod, 'Not MOD'
                    FA_SMILES_lst = []
                    FA_subchain_lst = []
                    FA_SMILES_str = ''

            if len(FA_subchain_lst) == 0:
                pass

            else:
                oxo_tail_lst = [')=O'] * len(FA_subchain_lst)
                FA_chain_lite_lst += oxo_tail_lst

            if FA_SMILES_lst != []:
                FA_SMILES_lst += FA_chain_lite_lst
                FA_SMILES_lst += FA_SMILES_end_lst
                FA_SMILES_str = ('').join(FA_SMILES_lst)

            else:
                FA_SMILES_str = ''

            print(FA_SMILES_str)

        if FA_SMILES_str == '':
            FA_SMILES_checker = 0

        else:
            FA_SMILES_checker = 1

        return FA_SMILES_str, FA_SMILES_checker, FA_ExactPos_checker


class GPabbr():

    def str2elem(self, hg=None, sn1=None, sn2=None, exact_pos=1):

        # set to use exact position of modifications
        if exact_pos == 1:
            usr_pos = 1

        if exact_pos == 0:
            usr_pos = 0

        else:
            usr_pos = 1

        mzcalc = Elem2Mass()
        # print self.usr_csv_pd.head(5)
        # print self.usr_csv_pd.columns.tolist()

        GP_elem_dct = {}

        # row = self.usr_csv_pd.iterrows()

        tmp_HG = hg
        tmp_R1 = sn1
        tmp_R2 = sn2

        # print tmp_HG, tmp_R1, tmp_R2
        HG_elem_dct = {}
        HG_SMILES_str = ''

        if tmp_HG == 'PA':
            HG_elem_dct['C'] = 3
            HG_elem_dct['H'] = 5
            HG_elem_dct['P'] = 1
            HG_elem_dct['O'] = 4
            HG_SMILES_str = 'OP(O)(OC[C@]([H])('

        if tmp_HG == 'PC':
            HG_elem_dct['C'] = 8
            HG_elem_dct['H'] = 16
            HG_elem_dct['N'] = 1
            HG_elem_dct['P'] = 1
            HG_elem_dct['O'] = 4
            HG_SMILES_str = '[O-]P(OCC[N+](C)(C)C)(OC[C@]([H])('

        if tmp_HG == 'PE':
            HG_elem_dct['C'] = 5
            HG_elem_dct['H'] = 10
            HG_elem_dct['N'] = 1
            HG_elem_dct['P'] = 1
            HG_elem_dct['O'] = 4
            HG_SMILES_str = 'OP(OCC[N])(OC[C@]([H])('

        if tmp_HG == 'PG':
            HG_elem_dct['C'] = 6
            HG_elem_dct['H'] = 11
            HG_elem_dct['P'] = 1
            HG_elem_dct['O'] = 6
            HG_SMILES_str = 'OP(OCC(O)CO)(OC[C@]([H])('

        if tmp_HG == 'PS':
            HG_elem_dct['C'] = 6
            HG_elem_dct['H'] = 10
            HG_elem_dct['N'] = 1
            HG_elem_dct['P'] = 1
            HG_elem_dct['O'] = 6
            HG_SMILES_str = 'OP(OCC(C(O)=O)N)(OC[C@]([H])('

        GP_elem_dct.update(HG_elem_dct)
        # print GP_elem_dct
        R1 = FAabbr(tmp_R1)
        (R1_elem_dct, R1_checker) = R1.FAstr2elem()
        (R1_SMILES_str, R1_SMILES_checker, R1_ExactPos_checker) = R1.FAstr2SMILES(exact_pos=usr_pos)

        R2 = FAabbr(tmp_R2)
        (R2_elem_dct, R2_checker) = R2.FAstr2elem()
        (R2_SMILES_str, R2_SMILES_checker, R2_ExactPos_checker) = R2.FAstr2SMILES(exact_pos=usr_pos)

        FA_checker = R1_checker + R2_checker

        FA_SMILES_checker = R1_SMILES_checker + R2_SMILES_checker

        FA_ExactPos_checker = R1_ExactPos_checker * R2_ExactPos_checker

        if FA_SMILES_checker == 0:
            GP_SMILES_str = ''

        elif FA_SMILES_checker == 1:
            GP_SMILES_str = ''

        elif FA_SMILES_checker == 2:
            if R1_SMILES_str != 'lyso' and R2_SMILES_str != 'lyso':
                GP_SMILES_lst = [HG_SMILES_str, R2_SMILES_str, ')C', R1_SMILES_str, ')=O']
                GP_SMILES_str = ('').join(GP_SMILES_lst)

            if R1_SMILES_str == 'lyso' and R2_SMILES_str != 'lyso':
                GP_SMILES_lst = [HG_SMILES_str, R2_SMILES_str, ')CO)=O']
                GP_SMILES_str = ('').join(GP_SMILES_lst)

            if R2_SMILES_str == 'lyso' and R1_SMILES_str != 'lyso':
                GP_SMILES_lst = [HG_SMILES_str, 'O)C', R1_SMILES_str, ')=O']
                GP_SMILES_str = ('').join(GP_SMILES_lst)

        if FA_checker == 2:
            for tmp_elem_k, tmp_elem_num in R1_elem_dct.items():
                # print tmp_mod_typ, tmp_elem_k, tmp_elem_num
                if tmp_elem_k in GP_elem_dct.keys():
                    GP_elem_dct[tmp_elem_k] += tmp_elem_num

                else:
                    GP_elem_dct[tmp_elem_k] = tmp_elem_num

            for tmp_elem_k, tmp_elem_num in R2_elem_dct.items():
                # print tmp_mod_typ, tmp_elem_k, tmp_elem_num
                if tmp_elem_k in GP_elem_dct.keys():
                    GP_elem_dct[tmp_elem_k] += tmp_elem_num

                else:
                    GP_elem_dct[tmp_elem_k] = tmp_elem_num

            GP_elem_str = ''
            GP_elem_idx_lst = ['C', 'H', 'O', 'P', 'N']

            for tmp_elem in GP_elem_idx_lst:
                if tmp_elem in GP_elem_dct.keys():
                    tmp_elem_num = GP_elem_dct[tmp_elem]
                    tmp_info = tmp_elem + str(tmp_elem_num)
                    GP_elem_str += tmp_info

                else:
                    pass

            # if tmp_R1 == '16_0' and tmp_R2 == '9_0 [1oxo]':
            # f1 = mzcalc.get_elem(GP_elem_str)
            #     mass = mzcalc.get_mass(f1)
            #     f = mzcalc.get_elem(GP_elem_str + 'H')
            #     mz = mzcalc.get_mass(f)
            #     print GP_elem_dct, GP_elem_str, mass, mz

            formula_org = mzcalc.get_elem(GP_elem_str)
            exact_mass = mzcalc.get_mass(formula_org)
            formula_H = mzcalc.get_elem(GP_elem_str + 'H')
            mz_H = mzcalc.get_mass(formula_H)
            formula_Na = mzcalc.get_elem(GP_elem_str + 'Na')
            mz_Na = mzcalc.get_mass(formula_Na)

        else:
            formula_org = ''
            exact_mass = ''
            mz_H = ''
            mz_Na = ''

        tmp_info_lst = [tmp_HG, tmp_R1, tmp_R2, formula_org, exact_mass, mz_H, mz_Na, GP_SMILES_str]
        print('SMILES Gnerated!')

        return GP_SMILES_str