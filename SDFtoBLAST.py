# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de
#

import re
import pandas as pd
from rdkit import Chem

from lpplibs.ExactMassCalc import Elem2Mass
from lpplibs.FormulaCalc import FormulaCalc

class AbbrFrag(object):

    def __init__(self, abbr, lpp_type):
        self.lpp_abbr = abbr
        self.lpp_type = lpp_type

        self.ox_dct = {}
        self.ox_dct['CH3'] = {'C': 1, 'H': 3}
        self.ox_dct['CH3COOH'] = {'C': 2, 'H': 4, 'O': 2}
        self.ox_dct['keto'] = {'H': -2, 'O': 1}
        self.ox_dct['C=O'] = {'H': -2, 'O': 1}
        self.ox_dct['OH'] = {'O': 1}
        self.ox_dct['CHO'] = {'H': -2, 'O': 1}
        self.ox_dct['COOH'] = {'H': -2, 'O': 2}
        # self.ox_dct['epoxy'] = {'H': -2, 'O': 1}
        # self.ox_dct['OOH'] = {'O': 2}

        # for PC
        self.pc_dct = {'C': 5, 'H': 14, 'O': 4, 'P': 1, 'N': 1}
        self.pc_gly_dct = {'C': 8, 'H': 20, 'O': 6, 'P': 1, 'N': 1}
        self.pc_gly_lyso_dct = {'C': 8, 'H': 18, 'O': 5, 'P': 1, 'N': 1}
        # for PE
        self.pe_dct = {'C': 2, 'H': 7, 'O': 4, 'P': 1, 'N': 1}
        self.pe_gly_lyso_dct = {'C': 5, 'H': 12, 'O': 5, 'P': 1, 'N': 1}
        # general
        self.glycerol_dct = {'C': 3, 'H': 8, 'O': 3}
        self.fa_dct = {'C': 1, 'H': 2, 'O': 2}
        self.water_dct = {'H': 2, 'O': 1}

    def get_mz(self, ion_df):

        mz_lst = []

        sn1_dct, sn2_dct, sn1_hydro_checker, sn2_hydro_checker = self.parse()

        mzcalc = Elem2Mass()
        formula_obj = FormulaCalc()
        # neutral exact mass/formula
        sn1_formula_txt = formula_obj.merge_dct('', sn1_dct)
        sn2_formula_txt = formula_obj.merge_dct('', sn2_dct)
        sn1_lyso_formula_txt = formula_obj.merge_dct(sn1_formula_txt, self.pc_gly_lyso_dct)
        sn2_lyso_formula_txt = formula_obj.merge_dct(sn2_formula_txt, self.pc_gly_lyso_dct)
        # sn1_lyso_formula_txt = formula_obj.merge_dct(sn1_formula_txt, self.pe_gly_lyso_dct)
        # sn2_lyso_formula_txt = formula_obj.merge_dct(sn2_formula_txt, self.pe_gly_lyso_dct)
        _m_txt = formula_obj.merge_dct(sn1_lyso_formula_txt, sn2_dct)
        # M
        m_formula_txt = formula_obj.merge_dct(_m_txt, {'H': -2, 'O': -1})
        m_mz = mzcalc.get_mass(m_formula_txt)
        # [M-H]-
        m_h_formula_txt = formula_obj.merge_dct(m_formula_txt, {'H': -1})
        m_h_mz = mzcalc.get_mass(m_h_formula_txt)
        # [M+FA-H]-
        m_fa_formula_txt = formula_obj.merge_dct(m_formula_txt, {'C': 1, 'H': 1, 'O': 2})
        m_fa_mz = mzcalc.get_mass(m_fa_formula_txt)
        print m_formula_txt, m_mz
        _pr_dct = {'M': (m_formula_txt, m_mz),
                  '[M-H]-': (m_h_formula_txt, m_h_mz),
                  '[M+FA-H]-': (m_fa_formula_txt, m_fa_mz)}

        for _ion_idx in ion_df.index.tolist():

            if self.lpp_type in ['OCP-alde', 'OAP']:
                pr_formula_txt = m_fa_formula_txt
                # pr_formula_txt = m_h_formula_txt

            elif self.lpp_type in ['OCP-acid']:
                pr_formula_txt = m_h_formula_txt

            if ion_df.loc[_ion_idx]['type'] == 'FRAG':
                if ion_df.loc[_ion_idx]['ion_name'] == 'C5H11NO5P-':
                    frag_mz = 196.0375
                    _frag_info = '%f %i \"%s\"\n' % (frag_mz, ion_df.loc[_ion_idx]['ion_i'], 'C5H11NO5P-')
                    mz_lst.append(_frag_info)
                if ion_df.loc[_ion_idx]['ion_name'] == '[sn1-H]-':
                    sn1_frag_formula_txt = formula_obj.merge_dct(sn1_formula_txt, {'H': -1})
                    sn1_mz = mzcalc.get_mass(sn1_frag_formula_txt)
                    print sn1_formula_txt, sn1_mz
                    _frag_info = '%f %i \"%s\"\n' % (sn1_mz, ion_df.loc[_ion_idx]['ion_i'], '[sn1-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[sn2-H]-':
                    sn2_frag_formula_txt = formula_obj.merge_dct(sn2_formula_txt, {'H': -1})
                    sn2_mz = mzcalc.get_mass(sn2_frag_formula_txt)
                    print sn2_formula_txt, sn2_mz
                    _frag_info = '%f %i \"%s\"\n' % (sn2_mz, ion_df.loc[_ion_idx]['ion_i'], '[sn2-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[sn2-H2O-H]-' and sn2_hydro_checker > 0:
                    _sn2_w_checker = 0
                    if self.lpp_type in ['OCP-alde', 'OCP-acid']:
                        _sn2_w_checker = 1
                    if sn2_hydro_checker == 0 and self.lpp_type == 'OAP':
                        _sn2_w_checker = 0
                    if sn2_hydro_checker > 0 and self.lpp_type == 'OAP':
                        _sn2_w_checker = 1
                    else:
                        pass

                    if _sn2_w_checker == 1:
                        sn2_frag_formula_txt = formula_obj.merge_dct(sn2_formula_txt, {'H': -3, 'O': -1})
                        sn2_mz = mzcalc.get_mass(sn2_frag_formula_txt)
                        print sn2_formula_txt, sn2_mz
                        _frag_info = '%f %i \"%s\"\n' % (sn2_mz, ion_df.loc[_ion_idx]['ion_i'], '[sn2-H2O-H]-')
                        mz_lst.append(_frag_info)

                # for OCP -COOH only
                if ion_df.loc[_ion_idx]['ion_name'] == '[sn2+CH2-H]-':

                    sn2_frag_formula_txt = formula_obj.merge_dct(sn2_formula_txt, {'H': 1, 'C': 1})
                    sn2_mz = mzcalc.get_mass(sn2_frag_formula_txt)
                    print sn2_formula_txt, sn2_mz
                    _frag_info = '%f %i \"%s\"\n' % (sn2_mz, ion_df.loc[_ion_idx]['ion_i'], '[sn2+CH2-H]-')
                    mz_lst.append(_frag_info)
                # for OCP -COOH only
                if ion_df.loc[_ion_idx]['ion_name'] == '[sn2-CH2O3-H]-':

                    sn2_frag_formula_txt = formula_obj.merge_dct(sn2_formula_txt, {'H': -3, 'C': -1, 'O': -3})
                    sn2_mz = mzcalc.get_mass(sn2_frag_formula_txt)
                    print sn2_formula_txt, sn2_mz
                    _frag_info = '%f %i \"%s\"\n' % (sn2_mz, ion_df.loc[_ion_idx]['ion_i'], '[sn2-CH2O3-H]-')
                    mz_lst.append(_frag_info)

            if ion_df.loc[_ion_idx]['type'] == 'NL':
                if ion_df.loc[_ion_idx]['ion_name'] == '[M-H]-' and self.lpp_type in ['OCP-alde', 'OCP-acid', 'OAP']:
                    _frag_formula_txt = formula_obj.merge_dct(pr_formula_txt, {'H': 0})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[M-CH3COOH-H]-' and self.lpp_type in ['OCP-alde', 'OAP']:
                    _frag_formula_txt = formula_obj.merge_dct(pr_formula_txt, {'C': -2, 'H': -4, 'O': -2})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-CH3COOH-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[M-PCsHG-H]-' and self.lpp_type in ['OCP-acid']:
                    _frag_formula_txt = formula_obj.merge_dct(pr_formula_txt, {'C': -3, 'H': -9, 'N': -1})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-PCsHG-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[M-sn2+H2O-H]-':
                    _frag_formula_txt = formula_obj.merge_dct(sn1_lyso_formula_txt, {'C': 0, 'H': -1})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-sn2+H2O-H]-')
                    mz_lst.append(_frag_info)
                if ion_df.loc[_ion_idx]['ion_name'] == '[M-sn2-CH3+H2O-H]-':
                    _frag_formula_txt = formula_obj.merge_dct(sn1_lyso_formula_txt, {'C': -1, 'H': -3})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-sn2-CH3+H2O-H]-')
                    mz_lst.append(_frag_info)
                if ion_df.loc[_ion_idx]['ion_name'] == '[M-sn2-CH3-H]-':
                    _frag_formula_txt = formula_obj.merge_dct(sn1_lyso_formula_txt, {'C': -1, 'H': -5, 'O': -1})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-sn2-CH3-H]-')
                    mz_lst.append(_frag_info)

                if ion_df.loc[_ion_idx]['ion_name'] == '[M-sn1-CH3+H2O-H]-':
                    _frag_formula_txt = formula_obj.merge_dct(sn2_lyso_formula_txt, {'C': -1, 'H': -3})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-sn1-CH3+H2O-H]-')
                    mz_lst.append(_frag_info)
                if ion_df.loc[_ion_idx]['ion_name'] == '[M-sn1-CH3-H]-':
                    _frag_formula_txt = formula_obj.merge_dct(sn2_lyso_formula_txt, {'C': -1, 'H': -5, 'O': -1})
                    _mz = mzcalc.get_mass(_frag_formula_txt)
                    print _frag_formula_txt, _mz
                    _frag_info = '%f %i \"%s\"\n' % (_mz, ion_df.loc[_ion_idx]['ion_i'], '[M-sn1-CH3-H]-')
                    mz_lst.append(_frag_info)

        mz_lst.append('\n')

        return mz_lst, _pr_dct

    def parse(self):

        oxpl_re = re.compile(r'(oxPC)([(])(.*)([/])(.*)([)])')

        _match = oxpl_re.match(self.lpp_abbr)
        if _match:
            _sn_lst = _match.groups()
            sn_lst = [_sn_lst[0], _sn_lst[2], _sn_lst[4]]
            print '_sn_lst', _sn_lst

            sn1_dct, sn1_hydro_checker = self.parse_sn(_sn_lst[2])
            try:
                sn2_dct, sn2_hydro_checker = self.parse_sn(_sn_lst[4])
            except:
                sn2_dct, sn2_hydro_checker = self.parse_sn(_sn_lst[2])

            print sn1_dct, sn2_dct, sn1_hydro_checker, sn2_hydro_checker

            return sn1_dct, sn2_dct, sn1_hydro_checker, sn2_hydro_checker

    def parse_sn(self, sn_abbr):

            sn_elem_dct = {'C': 0, 'H': 0, 'O': 2}

            # detect CHO or COOH
            alde_re = re.compile(r'(.*)(\(CHO@C\d{1,2}\))(.*)')
            acid_re = re.compile(r'(.*)(\(COOH@C\d{1,2}\))(.*)')
            alde_match = alde_re.match(sn_abbr)
            acid_match = acid_re.match(sn_abbr)
            if alde_match:
                print 'CHO'
                sn_elem_dct = {'C': 0, 'H': -4, 'O': 3}
                _sn_abbr_lst = alde_match.groups()
                _sn_abbr = ''.join([_sn_abbr_lst[0], _sn_abbr_lst[2]])
            elif acid_match:
                print 'COOH'
                sn_elem_dct = {'C': 0, 'H': -4, 'O': 4}
                _sn_abbr_lst = acid_match.groups()
                _sn_abbr = ''.join([_sn_abbr_lst[0], _sn_abbr_lst[2]])
            else:
                sn_elem_dct = {'C': 0, 'H': -2, 'O': 2}
                _sn_abbr = sn_abbr

            print '_sn_abbr', _sn_abbr, sn_abbr

            sn_unox_re = re.compile(r'(\d{1,2})([:])(\d)')
            sn_ox_re = re.compile(r'(\d{1,2})([:])(\d)([\[\(])(.*)([\)\]])')
            ox_re = re.compile(r'(\d)(x)(.{0,10})')

            _match_sn_unox = sn_unox_re.match(_sn_abbr)
            _match_sn_ox = sn_ox_re.match(_sn_abbr)

            _hydro_checker = 0

            if _match_sn_unox and not _match_sn_ox:
                _sn_lst = _match_sn_unox.groups()
                _re_sn = ''.join(_sn_lst)
                print '_sn_lst_unox', _sn_lst, _re_sn
                if _re_sn == _sn_abbr:
                    sn_elem_dct['C'] = int(_sn_lst[0])
                    sn_elem_dct['H'] += int(_sn_lst[0]) * 2 + 2 - 2 * int(_sn_lst[2])
                    print sn_abbr, _sn_lst, sn_elem_dct
                    return sn_elem_dct, _hydro_checker

            if _match_sn_ox:
                _sn_lst = _match_sn_ox.groups()
                print '_sn_lst_ox', _sn_lst, sn_elem_dct
                sn_elem_dct['C'] = int(_sn_lst[0])
                sn_elem_dct['H'] += int(_sn_lst[0]) * 2 + 2 - 2 * int(_sn_lst[2])
                _ox_info = _sn_lst[4]
                _ox_info_lst = _ox_info.split(',')
                print '_ox_info_lst', _ox_info_lst
                for _ox in _ox_info_lst:
                    _match_ox = ox_re.match(_ox)
                    if _match_ox:
                        _ox_lst = _match_ox.groups()
                        print '_ox_lst', _ox_lst
                        if _ox_lst[2] in self.ox_dct.keys():
                            if _ox_lst[2] == 'OH':
                                _hydro_checker += 1
                            _ox_count = int(_ox_lst[0])
                            _ox_dct = self.ox_dct[_ox_lst[2]]
                            print '_ox_dct', _ox_dct
                            for _elem in ['C', 'H', 'O']:
                                if _elem in _ox_dct.keys():
                                    _ox_add = _ox_count * _ox_dct[_elem]
                                    sn_elem_dct[_elem] += _ox_add
                return sn_elem_dct, _hydro_checker

sdf_file = 'oxPAPC_OH_KETO_sameDB.sdf'
msp_file = 'oxPAPC_OH_KETO_sameDB_2104.msp'
f_obj = file(msp_file, 'w')
f_obj.write('')
f_obj.close()

suppl = Chem.SDMolSupplier(sdf_file)

alde_re = re.compile(r'.*CHO.*')
acid_re = re.compile(r'.*COOH.*')
keto_re = re.compile(r'.*xketo.*')
hydro_re = re.compile(r'.*xOH.*')

for m in suppl:
    _id = m.GetProp('LM_ID')
    # _exactmass = m.GetProp('Exact_Mass')
    # _formula = m.GetProp('Formula')
    # _smiles = m.GetProp('SMILES')

    # give charge types [M+FA-H]- to OAP and OCP with -CHO, give [M-H]- to OCP with -COOH
    _charge_lst = ['[M-H]-', '[M+FA-H]-']

    if alde_re.match(_id):
        _charge_lst = ['[M+FA-H]-']
        # _charge_lst = ['[M-H]-']
        _frag_pattern = 'ion_scores_PLPC_OCP_CHO_neg.csv'
        _lpp_type = 'OCP-alde'
    elif acid_re.match(_id):
        _charge_lst = ['[M-H]-']
        _frag_pattern = 'ion_scores_PLPC_OCP_COOH_neg.csv'
        _lpp_type = 'OCP-acid'
    else:
        # _charge_lst = ['[M+FA-H]-']
        if hydro_re.match(_id) and keto_re.match(_id):
            _charge_lst = ['[M+FA-H]-']
            # _charge_lst = ['[M-H]-']
            _frag_pattern = 'ion_scores_PLPC_OAP_keto_OH_neg.csv'
            _lpp_type = 'OAP'
        elif hydro_re.match(_id) and not keto_re.match(_id):
            _charge_lst = ['[M+FA-H]-']
            # _charge_lst = ['[M-H]-']
            _frag_pattern = 'ion_scores_PLPC_OAP_OH_neg.csv'
            _lpp_type = 'OAP'
        elif keto_re.match(_id) and not hydro_re.match(_id):
            _charge_lst = ['[M+FA-H]-']
            # _charge_lst = ['[M-H]-']
            _frag_pattern = 'ion_scores_PLPC_OAP_keto_neg.csv'
            _lpp_type = 'OAP'
        else:
            _charge_lst = ['[M+FA-H]-']
            # _charge_lst = ['[M-H]-']
            _frag_pattern = 'ion_scores_PLPC_OAP_keto_neg.csv'
            _lpp_type = 'OAP'

    if _frag_pattern != '':

        # frag_obj = TheoFrag(_frag_pattern)
        _frag_df = pd.read_csv(_frag_pattern)

        for _charge in _charge_lst:

            # perform theo frag

            _theofrag = AbbrFrag(_id, _lpp_type)
            # try:
            mz_lst, pr_dct = _theofrag.get_mz(_frag_df)

            if _charge == '[M+FA-H]-':
                _exactmass = str(pr_dct['[M+FA-H]-'][1])
                _formula = pr_dct['M'][0]
            elif _charge == '[M-H]-':
                _exactmass = str(pr_dct['[M-H]-'][1])
                _formula = pr_dct['M'][0]

            else:
                _exactmass = 0.0
                _formula = 'N/A'

            _name = 'Name: %s (%s)\n' % (_id, _id)
            _id = 'LM_ID: %s\n' % _id
            _pr_type = 'Precursor_type: %s\n' % _charge
            _prmz = 'PrecursorMZ: %s\n' % _exactmass
            _formula = 'Formula: %s\n' % _formula
            # _comment = ('Comment: Parent=%s Mz_exact=%s ; %s; %s; %s; %s;\n' %
            #             (_exactmass, _exactmass, _id, _charge, _id, _formula))
            # _peaks_count = len(_frag_df.index.tolist())
            _peaks_count = len(mz_lst) - 1

            _num_peaks = 'Num Peaks: %i\n' % _peaks_count
            with open(msp_file, 'a') as msp_obj:
                _info_lst = [_name, _id, _pr_type, _prmz, _formula, _num_peaks]
                _info_lst.extend(mz_lst)
                # _info_lst.append('\n')
                msp_obj.writelines(_info_lst)


            # except TypeError:
            #     print 'skip!!!--->', _id
            #     pass

            break
