# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

from __future__ import division
import re


class Elem2Mass():
    # expands ((((M)N)O)P)Q to M*N*O*P*Q

    def get_elem(self, formula):
        self.formula = formula.encode("utf-8")
        while len(re.findall('\(\w*\)', str(self.formula))) > 0:
            parenthetical = re.findall('\(\w*\)[0-9]+', self.formula)
            for i in parenthetical:
                p = re.findall('[0-9]+', str(re.findall('\)[0-9]+', i)))
                j = re.findall('[A-Z][a-z]*[0-9]*', i)
                for n in range(0, len(j)):
                    numero = re.findall('[0-9]+', j[n])
                    if len(numero) != 0:
                        for k in numero:
                            nu = re.sub(k, str(int(int(k) * int(p[0]))), j[n])
                    else:
                        nu = re.sub(j[n], j[n] + p[0], j[n])
                    j[n] = nu
                newphrase = ''
                for m in j:
                    newphrase += str(m)
                self.formula = self.formula.replace(i, newphrase)
            if (len((re.findall('\(\w*\)[0-9]+', self.formula))) == 0) \
                    and (len(re.findall('\(\w*\)', self.formula)) != 0):
                self.formula = self.formula.replace('(', '')
                self.formula = self.formula.replace(')', '')

            print (self.formula)

        return self.formula

    #####################################################
    # Iterate over expanded formula to collect property
    #####################################################
    def get_mass(self, Formula_parsed):
        self.PeriodicTable = {
            'H': [1, 1, [1.0078250321, 2.0141017780], [0.999885, 0.0001157]],  # iupac '97 in water
            'D': [1, 1, [2.0141017780], [0.0001157]],  # iupac '97 in water
            'C': [6, 4, [12.0, 13.0033548378], [0.9893, 0.0107]],  # iupac '97
            'N': [7, 5, [14.0030740052, 15.0001088984], [0.99632, 0.00368]],  # iupac '97
            'O': [8, -2, [15.9949146221, 16.99913150, 17.9991604], [0.99757, 0.00038, 0.00205]],  # iupac '97
            'Na': [11, 1, [22.98976967], [1.0]],  # iupac '97
            'P': [15, 5, [30.97376151], [1.0]],  # iupac '97
            'S': [16, -2, [31.97207069, 32.97145850, 33.96786683, 35.96708088], [0.9493, 0.0076, 0.0429, 0.0002]],
            # iupac '97
            'Cl': [17, -1, [34.96885271, 36.96590260], [0.7578, 0.2422]],  # iupac '97
            'K': [19, 1, [38.9637069, 39.96399867, 40.96182597], [0.932581, 0.000117, 0.067302]],  # iupac '97
            'Ca': [20, 2, [39.9625912, 41.9586183, 42.9587668, 43.9554811, 45.9536928, 47.952534],
                   [0.96941, 0.00647, 0.00135, 0.02086, 0.00004, 0.00187]],  # iupac '97
            'Fe': [26, 3, [53.9396148, 55.9349421, 56.9353987, 57.9332805], [0.05845, 0.91754, 0.02119, 0.00282]],
            # iupac '97
            'Cu': [29, 2, [62.9296011, 64.9277937], [0.6917, 0.3083]],  # iupac '97
        }

        # print Formula_parsed
        ExactMass = 0
        while len(Formula_parsed) > 0:
            segments = re.findall('[A-Z][a-z]*[0-9]*', str(Formula_parsed))
            for i in range(0, len(segments)):
                # print 'i',i,segments[i]
                atom = re.findall('[A-Z][a-z]*', segments[i])
                number = re.findall('[0-9]+', segments[i])
                if len(number) == 0:
                    multiplier = 1
                else:
                    multiplier = float(number[0])
                # print self.PeriodicTable
                atomic_mass = self.PeriodicTable[atom[0]][2][0]
                seg_mass = atomic_mass * multiplier
                ExactMass += seg_mass
            Formula_parsed = re.sub(str(Formula_parsed), '', str(Formula_parsed))
        ExactMass = round(ExactMass, 5)
        return ExactMass

    def get_DBE(self, Formula_parsed):

        self.DBE_Table = {'H': -0.5, 'D': -0.5, 'C': 1, 'N': 0.5, 'O': 0, 'Na': -0.5, 'P': 0.5, 'S': 0,
                          'Cl': -0.5, 'K': -0.5, 'Ca': -1, 'Fe': 0, 'Cu': 0}
        DBE = 0.0
        while len(Formula_parsed) > 0:
            segments = re.findall('[A-Z][a-z]*[0-9]*', Formula_parsed)
            for i in range(0, len(segments)):
                # print 'i',i,segments[i]
                atom = re.findall('[A-Z][a-z]*', segments[i])
                number = re.findall('[0-9]+', segments[i])
                if len(number) == 0:
                    multiplier = 1
                else:
                    multiplier = float(number[0])
                atomic_DB = self.DBE_Table[atom[0]]
                seg_DB = atomic_DB * multiplier
                DBE += seg_DB
            Formula_parsed = re.sub(Formula_parsed, '', Formula_parsed)

        DBE += 1

        return DBE


class Elem2MZ(Elem2Mass):
    # ###############################################################################
    # expands ((((M)N)O)P)Q to M*N*O*P*Q
    # ###############################################################################

    def get_MZ_lst(self, formula, CHARGE):

        charge_lst = CHARGE.split(';')
        mz_lst = []
        calc = Elem2Mass()
        formula = calc.get_elem(formula)
        for C in charge_lst:
            if C.lower() == '1pos':
                temp_mz = calc.get_mass(formula)
                temp_mz += 1.0078250321
                temp_mz = round(temp_mz, 5)
                mz_lst.append(temp_mz)
            if C.lower() == '1neg':
                temp_mz = calc.get_mass(formula)
                temp_mz -= 1.0078250321
                temp_mz = round(temp_mz, 5)
                mz_lst.append(temp_mz)
            if C.lower() == '1na':
                temp_mz = calc.get_mass(formula)
                temp_mz += 22.98976967
                temp_mz = round(temp_mz, 5)
                mz_lst.append(temp_mz)
            if C.lower() == '1fa':
                temp_mz = calc.get_mass(formula)
                temp_mz += 1.0078250321 + 12.0 + 2 * 15.9949146221
                temp_mz = round(temp_mz, 5)
                mz_lst.append(temp_mz)
            if C.lower() == '2pos':
                temp_mz = calc.get_mass(formula)
                temp_mz = (temp_mz + 2 * 1.0078250321) / 2
                temp_mz = round(temp_mz, 5)
                mz_lst.append(temp_mz)

        return mz_lst


class MZcalc(object):
    def __init__(self):
        # iupac '97
        self.periodic_table_dct = {'H': [(1.0078250321, 0.999885), (2.0141017780, 0.0001157)],
                                   'D': [(2.0141017780, 0.0001157)],
                                   'C': [(12.0, 0.9893), (13.0033548378, 0.0107)],
                                   'N': [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
                                   'O': [(15.9949146221, 0.99757), (16.99913150, 0.00038), (17.9991604, 0.00205)],
                                   'Na': [(22.98976967, 1.0)],
                                   'P': [(30.97376151, 1.0)],
                                   'S': [(31.97207069, 0.9493), (32.97145850, 0.0076),
                                         (33.96786683, 0.0429), (35.96708088, 0.0002)],
                                   'K': [(38.9637069, 0.932581), (39.96399867, 0.000117), (40.96182597, 0.067302)],
                                   }

    def get_elements(self, formula):
        elem_dct = {}
        elem_key_lst = self.periodic_table_dct.keys()
        tmp_formula = formula

        elem_lst = re.findall('[A-Z][a-z]*[0-9]*', formula)
        for _e in range(0, len(elem_lst)):

            _elem = re.findall('[A-Z][a-z]*', elem_lst[_e])
            _elem_count = re.findall('[0-9]+', elem_lst[_e])
            if len(_elem_count) == 0:
                _elem_count = 1
            else:
                _elem_count = sum([int(x) for x in _elem_count])
            if _elem[0] in elem_dct.keys():
                elem_dct[_elem[0]] += _elem_count
            else:
                elem_dct[_elem[0]] = _elem_count

        return elem_dct

    def get_mono_mz(self, neutral_formula, charge):

        elem_dct = self.get_elements(neutral_formula)

        if charge in ['H-', '[M-H]-']:
            elem_dct['H'] += -1
        elif charge in ['+FA', '[M+FA-H]-', '[M+HCOO]-']:
            elem_dct['H'] += 1
            elem_dct['C'] += 1
            elem_dct['O'] += 2
        else:
            pass

        mono_mz = 0.0
        for _elem in elem_dct.keys():
            mono_mz += elem_dct[_elem] * self.periodic_table_dct[_elem][0][0]

        return mono_mz

# molecule='C42H72D9O8P1N1'
# CHG_lst = '1pos;1neg;2POS'
# #
# check1 = Elem2Mass()
# f = check1.get_elem(molecule)
# print f
# # DBE= check1.get_DBE(f)
# #
# # print DBE
# m1 = check1.get_mass(f)
# print 'm1',m1
#
# check = Elem2MZ()
# m = check.get_MZ_lst(molecule,CHG_lst)
# print m
#
# s = ''
# for x in m:
#    s = s + str(x)+';'
# s = s.strip(';')
