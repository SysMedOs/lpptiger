# -*- coding: utf-8 -*-
#
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPtiger.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#
from __future__ import division

try:  # python3
    import configparser
except NameError:  # python2
    import ConfigParser as configparser

import os
import re
import time

import xlrd
from PySide import QtCore, QtGui

from LPPtiger_UI import Ui_MainWindow
from LibHunter.Hunter_Core import huntlipids
from LibTheoLPP.TheoLPP_Core import theolpp


class LPPtiger_Main(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=os.getcwd()):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # current folder:
        if cwd is not None:
            print('User LPPtiger folder', cwd)
            self.theolpp_cwd = cwd
        else:
            auto_cwd = os.getcwd()
            print('User LPPtiger folder', auto_cwd)
            self.theolpp_cwd = auto_cwd

        # set groups for the radio buttons
        ox_level_group = QtGui.QButtonGroup(self.ui.tabWidget)
        ox_level_group.addButton(self.ui.ox_level1_rb)
        ox_level_group.addButton(self.ui.ox_level2_rb)
        ox_level_group.addButton(self.ui.ox_level3_rb)

        prostane_group = QtGui.QButtonGroup(self.ui.tabWidget)
        prostane_group.addButton(self.ui.prostane_yes_rb)
        prostane_group.addButton(self.ui.prostane_no_rb)

        prostane_ox_group = QtGui.QButtonGroup(self.ui.tabWidget)
        prostane_ox_group.addButton(self.ui.prostane_ox_yes_rb)
        prostane_ox_group.addButton(self.ui.prostane_ox_no_rb)

        msp_group = QtGui.QButtonGroup(self.ui.tabWidget)
        msp_group.addButton(self.ui.spectra_yes_rb)
        msp_group.addButton(self.ui.spectra_no_rb)

        # slots for tab a
        # hide the msp part
        self.ui.label_14.hide()
        self.ui.exceltab_cb.hide()
        # self.ui.label_12.hide()
        # self.ui.prostane_ox_yes_rb.hide()
        # self.ui.prostane_ox_no_rb.hide()
        QtCore.QObject.connect(self.ui.spectra_yes_rb, QtCore.SIGNAL("clicked()"), self.set_spec_t)
        QtCore.QObject.connect(self.ui.spectra_no_rb, QtCore.SIGNAL("clicked()"), self.set_spec_f)
        QtCore.QObject.connect(self.ui.prostane_yes_rb, QtCore.SIGNAL("clicked()"), self.set_prostane_t)
        QtCore.QObject.connect(self.ui.prostane_no_rb, QtCore.SIGNAL("clicked()"), self.set_prostane_f)
        QtCore.QObject.connect(self.ui.load_lipid_pb, QtCore.SIGNAL("clicked()"), self.load_lipid_list)
        QtCore.QObject.connect(self.ui.save_sdf_pb, QtCore.SIGNAL("clicked()"), self.save_sdf)
        QtCore.QObject.connect(self.ui.save_msp_pb, QtCore.SIGNAL("clicked()"), self.save_msp)
        QtCore.QObject.connect(self.ui.run_theolpp_pb, QtCore.SIGNAL("clicked()"), self.run_theolpp)
        QtCore.QObject.connect(self.ui.max_ox_spb, QtCore.SIGNAL("valueChanged(int)"), self.set_ox_max)

        # slots for tab b
        QtCore.QObject.connect(self.ui.tab_b_loadlpppath_pb, QtCore.SIGNAL("clicked()"), self.b_load_sum_sdf)
        QtCore.QObject.connect(self.ui.tab_b_loadfapath_pb, QtCore.SIGNAL("clicked()"), self.b_load_sum_fa)
        # QtCore.QObject.connect(self.ui.tab_b_loadsdfpath_pb, QtCore.SIGNAL("clicked()"), self.b_load_sdf)
        # QtCore.QObject.connect(self.ui.tab_b_loadmsppath_pb, QtCore.SIGNAL("clicked()"), self.b_load_msp)
        QtCore.QObject.connect(self.ui.tab_b_ms2mzml_pb, QtCore.SIGNAL("clicked()"), self.b_load_mzml)
        QtCore.QObject.connect(self.ui.tab_b_saveimgfolder_pb, QtCore.SIGNAL("clicked()"), self.b_save_img2folder)
        QtCore.QObject.connect(self.ui.tab_b_sumxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.b_save_output)
        QtCore.QObject.connect(self.ui.tab_b_runhunter_pb, QtCore.SIGNAL("clicked()"), self.b_run_hunter)

        # slots for tab c
        QtCore.QObject.connect(self.ui.tab_c_mod_lst_pb, QtCore.SIGNAL("clicked()"), self.set_general_mod)
        QtCore.QObject.connect(self.ui.tab_c_fa_lst_pb, QtCore.SIGNAL("clicked()"), self.set_fa_white_lst)
        QtCore.QObject.connect(self.ui.tab_c_prostane_mod_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_mod)
        QtCore.QObject.connect(self.ui.tab_c_prostane_abbr_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_abbr)
        QtCore.QObject.connect(self.ui.tab_c_frag_pattern_pb, QtCore.SIGNAL("clicked()"), self.set_frag_pattern)
        QtCore.QObject.connect(self.ui.tab_c_hgcfg_pb, QtCore.SIGNAL("clicked()"), self.set_hg_specifc_pattern)
        QtCore.QObject.connect(self.ui.tab_c_scorecfg_pb, QtCore.SIGNAL("clicked()"), self.set_wfrag)
        QtCore.QObject.connect(self.ui.set_default_pb, QtCore.SIGNAL("clicked()"), self.set_default_cfg)
        QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"), self.calc_snr_amp)

        # load configurations
        config = configparser.ConfigParser()
        config.read('configure.ini')
        if config.has_section('settings'):
            user_cfg = 'settings'
        else:
            if config.has_section('default'):
                user_cfg = 'default'
            else:
                user_cfg = ''
        if len(user_cfg) > 2:
            options = config.options(user_cfg)
            if 'general_mod_lst' in options:
                self.ui.tab_c_mod_lst_le.setText(config.get(user_cfg, 'general_mod_lst'))
            if 'fa_white_list' in options:
                self.ui.tab_c_fa_lst_le.setText(config.get(user_cfg, 'fa_white_list'))
            if 'prostane_mod_lst' in options:
                self.ui.tab_c_prostane_mod_lst_le.setText(config.get(user_cfg, 'prostane_mod_lst'))
            if 'prostane_abbr' in options:
                self.ui.tab_c_prostane_abbr_lst_le.setText(config.get(user_cfg, 'prostane_abbr'))
            if 'frag_patterns' in options:
                self.ui.tab_c_frag_pattern_le.setText(config.get(user_cfg, 'frag_patterns'))
            if 'hg_specifc_lst' in options:
                self.ui.tab_c_hgcfg_le.setText(config.get(user_cfg, 'hg_specifc_lst'))
            if 'wfrag_lst' in options:
                self.ui.tab_c_scorecfg_le.setText(config.get(user_cfg, 'wfrag_lst'))
            if 'wfrag_lst' in options:
                self.ui.tab_c_snratio_spb.setValue(int(config.get(user_cfg, 'sn_ratio')))

    def set_spec_t(self):
        self.ui.save_msp_le.show()
        self.ui.save_msp_pb.show()
        self.ui.save_msp_lb.show()

    def set_spec_f(self):
        self.ui.save_msp_le.hide()
        self.ui.save_msp_pb.hide()
        self.ui.save_msp_lb.hide()

    def set_prostane_t(self):
        self.ui.label_12.show()
        self.ui.prostane_ox_yes_rb.show()
        self.ui.prostane_ox_no_rb.show()

    def set_prostane_f(self):
        self.ui.label_12.hide()
        self.ui.prostane_ox_yes_rb.hide()
        self.ui.prostane_ox_no_rb.hide()

    def set_ox_max(self):
        _ox_max = self.ui.max_ox_spb.value()
        self.ui.max_keto_spb.setMaximum(_ox_max)

    def set_general_mod(self):
        self.ui.tab_c_mod_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_mod_lst_le.setText(a_load_table_str)

    def set_fa_white_lst(self):
        self.ui.tab_c_fa_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_fa_lst_le.setText(a_load_table_str)

    def set_prostane_mod(self):
        self.ui.tab_c_prostane_mod_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_prostane_mod_lst_le.setText(a_load_table_str)

    def set_prostane_abbr(self):
        self.ui.tab_c_prostane_abbr_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_prostane_abbr_lst_le.setText(a_load_table_str)

    def set_frag_pattern(self):
        self.ui.tab_c_frag_pattern_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_frag_pattern_le.setText(a_load_table_str)

    def set_hg_specifc_pattern(self):
        self.ui.tab_c_hgcfg_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_hgcfg_le.setText(a_load_table_str)

    def set_wfrag(self):
        self.ui.tab_c_scorecfg_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.tab_c_scorecfg_le.setText(a_load_table_str)

    def set_default_cfg(self):
        config = configparser.ConfigParser()
        with open('configure.ini', 'w') as default_cfg:
            config.add_section('settings')
            config.set('settings', 'general_mod_lst', self.ui.tab_c_mod_lst_le.text())
            config.set('settings', 'fa_white_list', self.ui.tab_c_fa_lst_le.text())
            config.set('settings', 'prostane_mod_lst', self.ui.tab_c_prostane_mod_lst_le.text())
            config.set('settings', 'prostane_abbr', self.ui.tab_c_prostane_abbr_lst_le.text())
            config.set('settings', 'frag_patterns', self.ui.tab_c_frag_pattern_le.text())
            config.set('settings', 'hg_specifc_lst', self.ui.tab_c_hgcfg_le.text())
            config.set('settings', 'wfrag_lst', self.ui.tab_c_scorecfg_le.text())
            config.set('settings', 'sn_ratio', str(self.ui.tab_c_snratio_spb.value()))
            config.write(default_cfg)

    def load_lipid_list(self):
        self.ui.load_lipid_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters(['Tables (*.xlsx)'])
        load_table_dialog.selectNameFilter('Tables (*.xlsx)')
        if load_table_dialog.exec_():
            self.ui.label_14.show()
            self.ui.exceltab_cb.show()
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.load_lipid_le.setText(a_load_table_str)
            lipid_xlsx = xlrd.open_workbook(a_load_table_str)

            sheetnames_lst = lipid_xlsx.sheet_names()
            self.ui.exceltab_cb.addItems(sheetnames_lst)

    def save_sdf(self):
        save_sdf_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.sdf')
        self.ui.save_sdf_le.clear()
        save_sdf_str = os.path.abspath(save_sdf_path[0])
        self.ui.save_sdf_le.setText(save_sdf_str)

    def save_msp(self):
        save_msp_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.msp')
        self.ui.save_msp_le.clear()
        save_msp_str = os.path.abspath(save_msp_path[0])
        self.ui.save_msp_le.setText(save_msp_str)

    def run_theolpp(self):

        self.ui.run_status_te.append('start!')
        ox_level = 1
        if self.ui.ox_level1_rb.isChecked():
            ox_level = 1
        elif self.ui.ox_level2_rb.isChecked():
            ox_level = 2
        elif self.ui.ox_level3_rb.isChecked():
            ox_level = 3

        oap_mode, ocp_mode, lyso_oap_mode, lyso_ocp_mode = 0, 0, 0, 0
        if self.ui.oap_chb.isChecked():
            oap_mode = 1
        if self.ui.ocp_chb.isChecked():
            ocp_mode = 1
        if self.ui.lyso_oap_chb.isChecked():
            lyso_oap_mode = 1
        if self.ui.lyso_ocp_chb.isChecked():
            lyso_ocp_mode = 1

        ox_max = self.ui.max_ox_spb.value()
        keto_max = self.ui.max_keto_spb.value()
        ooh_max = self.ui.max_ooh_spb.value()
        epoxy_max = self.ui.max_epoxy_spb.value()

        prostane_mode = 0
        if self.ui.prostane_yes_rb.isChecked():
            prostane_mode = 1
        elif self.ui.prostane_no_rb.isChecked():
            prostane_mode = 0

        ox_prostane_mode = 0
        if self.ui.prostane_ox_yes_rb.isChecked():
            ox_prostane_mode = 1
        elif self.ui.prostane_ox_no_rb.isChecked():
            ox_prostane_mode = 0

        msp_mode = 1  # by default save msp is checked
        if self.ui.spectra_yes_rb.isChecked():
            msp_mode = 1

        elif self.ui.spectra_no_rb.isChecked():
            msp_mode = 0

        lipid_lst_path = self.ui.load_lipid_le.text()

        lipid_tab = self.ui.exceltab_cb.currentText()
        lipid_class = self.ui.pl_class_cb.currentText()
        lipid_class = lipid_class[-3:-1]  # get PL abbr

        sdf_path = self.ui.save_sdf_le.text()
        msp_path = self.ui.save_msp_le.text()
        mod_lst_path = self.ui.tab_c_mod_lst_le.text()
        fa_lst_path = self.ui.tab_c_fa_lst_le.text()
        prostane_mod_path = self.ui.tab_c_prostane_mod_lst_le.text()
        prostane_abbr_path = self.ui.tab_c_prostane_abbr_lst_le.text()
        frag_pattern_path = self.ui.tab_c_frag_pattern_le.text()
        hgcfg_path = self.ui.tab_c_hgcfg_le.text()

        param_dct = {'lipid_class': lipid_class, 'ox_level': ox_level,
                     'oap_mode': oap_mode, 'ocp_mode': ocp_mode,
                     'lyso_oap_mode': lyso_oap_mode, 'lyso_ocp_mode': lyso_ocp_mode,
                     'ox_max': ox_max, 'keto_max': keto_max, 'ooh_max': ooh_max, 'epoxy_max': epoxy_max,
                     'lipid_lst_path': lipid_lst_path, 'lipid_tab': lipid_tab,
                     'prostane_mode': prostane_mode, 'ox_prostane_mode': ox_prostane_mode,
                     'sdf_path': sdf_path, 'msp_mode': msp_mode, 'msp_path': msp_path,
                     'mod_lst_path': mod_lst_path, 'fa_lst_path': fa_lst_path, 'prostane_mod_path': prostane_mod_path,
                     'prostane_abbr_path': prostane_abbr_path, 'frag_pattern_path': frag_pattern_path,
                     'pl_hg_path': hgcfg_path}

        print(param_dct)
        info_1, info_2 = theolpp(param_dct)
        self.ui.run_status_te.append(info_1)
        self.ui.run_status_te.append(info_2)
        self.ui.run_status_te.append('Finished!\n !Please exit LPPtiger before open the .sdf output file!')

        # try:
        #     info_1, info_2 = theolpp(param_dct)
        #     self.ui.run_status_te.append(info_1)
        #     self.ui.run_status_te.append(info_2)
        #     self.ui.run_status_te.append('Finished!\n !Please exit LPPtiger before open the .sdf output file!')
        #
        # except:
        #     self.ui.run_status_te.append('!! An error has occurred, please check your settings !!')

    def b_load_sum_sdf(self):
        b_load_lipidstable_dialog = QtGui.QFileDialog(self)
        b_load_lipidstable_dialog.setNameFilters(['MS Excel files (*.xlsx *.XLSX)'])
        b_load_lipidstable_dialog.selectNameFilter('MS Excel files (*.xlsx *.XLSX)')
        if b_load_lipidstable_dialog.exec_():
            self.ui.tab_b_loadlpppath_le.clear()
            b_load_xlsx_str = b_load_lipidstable_dialog.selectedFiles()[0]
            b_load_xlsx_str = os.path.abspath(b_load_xlsx_str)
            self.ui.tab_b_loadlpppath_le.setText(b_load_xlsx_str)

    def b_load_sum_fa(self):
        b_load_lipidstable_dialog = QtGui.QFileDialog(self)
        b_load_lipidstable_dialog.setNameFilters(['MS Excel files (*.xlsx *.XLSX)'])
        b_load_lipidstable_dialog.selectNameFilter('MS Excel files (*.xlsx *.XLSX)')
        if b_load_lipidstable_dialog.exec_():
            self.ui.tab_b_loadfapath_le.clear()
            b_load_xlsx_str = b_load_lipidstable_dialog.selectedFiles()[0]
            b_load_xlsx_str = os.path.abspath(b_load_xlsx_str)
            self.ui.tab_b_loadfapath_le.setText(b_load_xlsx_str)
    #
    # def b_load_sdf(self):
    #     b_load_lipidstable_dialog = QtGui.QFileDialog(self)
    #     b_load_lipidstable_dialog.setNameFilters([u'SDF files (*.sdf *.SDF)'])
    #     b_load_lipidstable_dialog.selectNameFilter(u'SDF files (*.sdf *.SDF)')
    #     if b_load_lipidstable_dialog.exec_():
    #         self.ui.tab_b_loadsdfpath_le.clear()
    #         b_load_sdf_str = b_load_lipidstable_dialog.selectedFiles()[0]
    #         b_load_sdf_str = os.path.abspath(b_load_sdf_str)
    #         self.ui.tab_b_loadsdfpath_le.setText(unicode(b_load_sdf_str))
    #
    # def b_load_msp(self):
    #     b_load_lipidstable_dialog = QtGui.QFileDialog(self)
    #     b_load_lipidstable_dialog.setNameFilters([u'SDF files (*.msp *.MSP)'])
    #     b_load_lipidstable_dialog.selectNameFilter(u'SDF files (*.msp *.MSP)')
    #     if b_load_lipidstable_dialog.exec_():
    #         self.ui.tab_b_loadmsppath_le.clear()
    #         b_load_msp_str = b_load_lipidstable_dialog.selectedFiles()[0]
    #         b_load_msp_str = os.path.abspath(b_load_msp_str)
    #         self.ui.tab_b_loadmsppath_le.setText(unicode(b_load_msp_str))

    def b_load_mzml(self):
        b_load_mzml_dialog = QtGui.QFileDialog(self)
        b_load_mzml_dialog.setNameFilters(['mzML spectra files (*.mzML *.mzml)'])
        b_load_mzml_dialog.selectNameFilter('mzML spectra files (*.mzML *.mzml)')
        if b_load_mzml_dialog.exec_():
            self.ui.tab_b_ms2mzml_le.clear()
            b_load_mzml_str = b_load_mzml_dialog.selectedFiles()[0]
            b_load_mzml_str = os.path.abspath(b_load_mzml_str)
            self.ui.tab_b_ms2mzml_le.setText(b_load_mzml_str)

    def b_save_img2folder(self):
        self.ui.tab_b_saveimgfolder_le.clear()
        e_save_img2folder_str = QtGui.QFileDialog.getExistingDirectory()
        self.ui.tab_b_saveimgfolder_le.setText(e_save_img2folder_str)

    def b_save_output(self):
        e_save_output_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.xlsx')
        self.ui.tab_b_sumxlsxpath_le.clear()
        e_save_output_str = os.path.abspath(e_save_output_path[0])
        self.ui.tab_b_sumxlsxpath_le.setText(e_save_output_str)

    def b_run_hunter(self):

        # self.ui.tab_b_loadlpppath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_lipids\Lv1_max3O_max1Keto_prostane\PC_Lv1_max3O_max1Keto_prostane.xlsx')
        # self.ui.tab_b_loadfapath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_lipids\Lv1_max3O_max1Keto_prostane\PC_Lv1_max3O_max1Keto_prostane.xlsx')
        # # self.ui.tab_b_loadlpppath_le.setText(r'D:\Project_lpptiger\output_sdf\PLstd\PC_std.xlsx')
        # # self.ui.tab_b_loadfapath_le.setText(r'D:\Project_lpptiger\output_sdf\PLstd\PC_std_FA_SUM.xlsx')
        # # self.ui.tab_b_loadsdfpath_le.setText(r'D:\Project_lpptiger\output_sdf\PC_FP.sdf')
        # # self.ui.tab_b_loadmsppath_le.setText(r'D:\Project_lpptiger\output_sdf\PC_FP.msp')
        # self.ui.tab_b_ms2mzml_le.setText(r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_I.mzML')
        # # self.ui.tab_b_ms2mzml_le.setText(r'D:\Project_lpptiger\mzML\131015_PLPC_400ng_new_neg_LMQ15.mzML')
        # # self.ui.tab_b_ms2mzml_le.setText(r'D:\Synapt_rawspectra\oxPLstd\180816_oxPC_10ng.mzML')
        # # self.ui.tab_b_saveimgfolder_le.setText(r'D:\Project_lpptiger\output_sdf\hunter_output\')
        # self.ui.tab_b_saveimgfolder_le.setText(r'D:\Project_lpptiger\output_sdf\CM_LPPs\PC')
        # # self.ui.tab_b_sumxlsxpath_le.setText(r'D:\Project_lpptiger\output_sdf\hunter_output\test_PC_FP_CM_C18.xlsx')
        # self.ui.tab_b_sumxlsxpath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_LPPs\PC\CM_PC_Lv1_prostane.xlsx')
        self.ui.tab_b_rtstart_dspb.setValue(3)
        self.ui.tab_b_rtend_dspb.setValue(31)
        self.ui.tab_b_msppm_spb.setValue(20)
        self.ui.tab_b_ms2ppm_spb.setValue(50)
        self.ui.tab_b_hgppm_spb.setValue(200)
        self.ui.tab_b_dda_spb.setValue(12)
        self.ui.tab_b_msthreshold_spb.setValue(1000)
        self.ui.tab_b_ms2threshold_spb.setValue(10)
        self.ui.tab_b_score_spb.setValue(20.5)
        self.ui.tab_b_isotopescore_spb.setValue(80.0)
        self.ui.tab_b_mzstart_dspb.setValue(600.0)
        self.ui.tab_b_mzend_dspb.setValue(1000.0)
        self.ui.tab_b_ms2infoth_dspb.setValue(1)
        self.ui.tab_b_ms2hginfoth_dspb.setValue(1)

        if self.ui.vendor_waters_rb.isChecked():
            usr_vendor = 'waters'
        elif self.ui.vendor_thermo_rb.isChecked():
            usr_vendor = 'thermo'
        else:
            usr_vendor = 'waters'
        if self.ui.mode_lcms_rb.isChecked():
            usr_exp_mode = 'LC-MS'
        elif self.ui.mode_static_rb.isChecked():
            usr_exp_mode = 'Static-MS'
        else:
            usr_exp_mode = 'LC-MS'
        print('Vendor mode = %s, Experiment mode = %s' % (usr_vendor, usr_exp_mode))
        print('Hunter started!')
        _pl_class_info = str(self.ui.tab_b_lipidclass_cmb.currentText())

        pl_class_checker = re.compile(r'(.*)( [\(])(\w{2,3})([\)] )(.*)')

        pl_class_match = pl_class_checker.match(_pl_class_info)

        if pl_class_match:
            pl_class_info_lst = pl_class_match.groups()
            _pl_class = pl_class_info_lst[2]
            _pl_charge = pl_class_info_lst[4]
        else:
            _pl_class = 'PC'
            _pl_charge = '[M+HCOO]-'

        lpp_sum_info_path_str = str(self.ui.tab_b_loadlpppath_le.text())
        fa_sum_path_str = str(self.ui.tab_b_loadfapath_le.text())
        # sdf_path_str = str(self.ui.tab_b_loadsdfpath_le.text())
        # msp_path_str = str(self.ui.tab_b_loadmsppath_le.text())
        mzml_path_str = str(self.ui.tab_b_ms2mzml_le.text())
        img_output_folder_str = str(self.ui.tab_b_saveimgfolder_le.text())
        xlsx_output_path_str = str(self.ui.tab_b_sumxlsxpath_le.text())

        rt_start = self.ui.tab_b_rtstart_dspb.value()
        rt_end = self.ui.tab_b_rtend_dspb.value()
        mz_start = self.ui.tab_b_mzstart_dspb.value()
        mz_end = self.ui.tab_b_mzend_dspb.value()
        dda_top = self.ui.tab_b_dda_spb.value()
        ms_th = self.ui.tab_b_msthreshold_spb.value()
        ms2_th = self.ui.tab_b_ms2threshold_spb.value()
        ms_ppm = self.ui.tab_b_msppm_spb.value()
        ms2_ppm = self.ui.tab_b_ms2ppm_spb.value()
        pr_window = self.ui.tab_b_prwindow_spb.value()
        hg_th = self.ui.tab_b_hgthreshold_spb.value()
        hg_ppm = self.ui.tab_b_hgppm_spb.value()
        score_filter = self.ui.tab_b_score_spb.value()
        isotope_score_filter = self.ui.tab_b_isotopescore_spb.value()
        ms2_info_threshold = self.ui.tab_b_ms2infoth_dspb.value() * 0.01
        hgms2_info_threshold = self.ui.tab_b_ms2hginfoth_dspb.value() * 0.01

        # from settings tab
        lipid_specific_cfg = self.ui.tab_c_hgcfg_le.text()
        score_cfg = self.ui.tab_c_scorecfg_le.text()
        sn_ratio = self.ui.tab_c_snratio_spb.value()

        try:
            if os.path.isfile(lipid_specific_cfg):
                pass
        except IOError:
            self.ui.tab_b_statusrun_pte.insertPlainText('!! Failed to load configuration for Phospholipids ',
                                                        'specific signal, please check your settings!!')
        try:
            if os.path.isfile(score_cfg):
                pass
        except IOError:
            self.ui.tab_b_statusrun_pte.insertPlainText('!! Failed to load score configuration',
                                                        'please check your settings!!')

        usr_score_mode = self.ui.tab_b_scoremode_cmb.currentIndex()
        if usr_score_mode == 0:
            print(self.ui.tab_b_scoremode_cmb.currentText())
            rank_score = True
        else:
            print(self.ui.tab_b_scoremode_cmb.currentText())
            rank_score = False

        usr_isotope_score_mode = self.ui.tab_b_isotopescoremode_cmb.currentIndex()
        if usr_isotope_score_mode == 0:
            print(self.ui.tab_b_isotopescoremode_cmb.currentText())
            fast_isotope = False
        else:
            print(self.ui.tab_b_isotopescoremode_cmb.currentText())
            fast_isotope = True

        start_time_str = time.strftime("%Y-%m-%d_%H-%M", time.localtime())

        # 'sdf_path_str': sdf_path_str,
        # 'msp_path_str': msp_path_str,

        hunter_param_dct = {'lpp_sum_info_path_str': lpp_sum_info_path_str,
                            'fa_sum_path_str': fa_sum_path_str,
                            'mzml_path_str': mzml_path_str,
                            'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str, 'rt_start': rt_start, 'rt_end': rt_end,
                            'mz_start': mz_start, 'mz_end': mz_end, 'dda_top': dda_top, 'ms_th': ms_th,
                            'ms2_th': ms2_th, 'ms_ppm': ms_ppm, 'ms2_ppm': ms2_ppm, 'hg_th': hg_th, 'hg_ppm': hg_ppm,
                            'score_filter': score_filter, 'isotope_score_filter': isotope_score_filter,
                            'lipid_type': _pl_class, 'charge_mode': _pl_charge, 'pr_window': pr_window,
                            'lipid_specific_cfg': lipid_specific_cfg, 'score_cfg': score_cfg, 'vendor': usr_vendor,
                            'ms2_infopeak_threshold': ms2_info_threshold,
                            'ms2_hginfopeak_threshold': hgms2_info_threshold,
                            'rank_score': rank_score, 'fast_isotope': fast_isotope, 'sn_ratio': sn_ratio,
                            'hunter_folder': self.theolpp_cwd,
                            'hunter_start_time': start_time_str, 'Experiment_mode': usr_exp_mode}

        param_log_output_path_str = (str(self.ui.tab_b_saveimgfolder_le.text()) +
                                     '/LPPtiger_Params-Log_%s.txt' % start_time_str
                                     )

        config = configparser.ConfigParser()
        with open(param_log_output_path_str, 'w') as usr_param_cfg:
            config.add_section('parameters')
            for param in hunter_param_dct.keys():
                config.set('parameters', param, str(hunter_param_dct[param]))
            config.write(usr_param_cfg)

        print(hunter_param_dct)

        tot_run_time = huntlipids(hunter_param_dct)

        if isinstance(tot_run_time, str):
            self.ui.tab_b_statusrun_pte.insertPlainText(tot_run_time)
            self.ui.tab_b_statusrun_pte.insertPlainText('>>> >>> >>> FINISHED <<< <<< <<<')

        else:
            self.ui.tab_b_statusrun_pte.insertPlainText('!! Failed read input files, please check vendor settings!!')

    def calc_snr_amp(self):
        import math

        usr_max_sn_ratio = self.ui.tab_c_testsnratio_spb.value()
        usr_test_sn_ratio = self.ui.tab_c_trysnratio_spb.value()

        usr_amp_factor = 100 / (20 * math.log10(usr_max_sn_ratio))

        origin_snr_score = 20 * math.log10(usr_test_sn_ratio)
        amp_snr_score = usr_amp_factor * 20 * math.log10(usr_test_sn_ratio)
        if amp_snr_score >= 100.0:
            amp_snr_score = 100.0

        self.ui.tab_c_originsnr_le.setText('%.1f' % origin_snr_score)
        self.ui.tab_c_ampsnr_le.setText('%.1f' % amp_snr_score)

if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = LPPtiger_Main()
    window.show()
    sys.exit(app.exec_())
