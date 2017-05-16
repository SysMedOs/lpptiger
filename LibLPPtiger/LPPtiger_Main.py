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

import ConfigParser as configparser

import os
import glob
import re
import time
import multiprocessing
import multiprocessing.pool

import xlrd
from PySide import QtCore, QtGui

from LPPtiger_UI import Ui_MainWindow
from LibHunter.Hunter_Core import huntlipids
from LibTheoLPP.TheoLPP_Core import theolpp


class LPPtigerMain(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, cwd=os.getcwd()):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.tabWidget.setCurrentIndex(0)
        self.ui.tabWidget.removeTab(4)

        self.ui.version_lb.setText('LPPtiger Beta Version: 16, May, 2017')

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

        msp_group = QtGui.QButtonGroup(self.ui.tabWidget)
        msp_group.addButton(self.ui.spectra_yes_rb)
        msp_group.addButton(self.ui.spectra_no_rb)

        self.ui.pl_class_cb.removeItem(8)
        self.ui.pl_class_cb.removeItem(7)
        self.ui.pl_class_cb.removeItem(6)
        self.ui.pl_class_cb.removeItem(4)

        self.ui.tab_b_lipidclass_cmb.removeItem(8)
        self.ui.tab_b_lipidclass_cmb.removeItem(7)
        self.ui.tab_b_lipidclass_cmb.removeItem(6)
        self.ui.tab_b_lipidclass_cmb.removeItem(4)

        self.d_set_multi_mode()

        # slots for tab a
        # hide the msp part
        self.ui.label_14.hide()
        self.ui.exceltab_cb.hide()
        self.ui.tab_b_ms1max_lb.hide()
        self.ui.tab_b_ms1max_spb.hide()
        # self.ui.label_12.hide()
        # self.ui.prostane_ox_yes_rb.hide()
        # self.ui.prostane_ox_no_rb.hide()
        QtCore.QObject.connect(self.ui.spectra_yes_rb, QtCore.SIGNAL("clicked()"), self.set_spec_t)
        QtCore.QObject.connect(self.ui.spectra_no_rb, QtCore.SIGNAL("clicked()"), self.set_spec_f)
        # QtCore.QObject.connect(self.ui.prostane_yes_rb, QtCore.SIGNAL("clicked()"), self.set_prostane_t)
        # QtCore.QObject.connect(self.ui.prostane_no_rb, QtCore.SIGNAL("clicked()"), self.set_prostane_f)
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
        QtCore.QObject.connect(self.ui.tab_b_ms1max_chb, QtCore.SIGNAL("clicked()"), self.set_ms1_max)
        QtCore.QObject.connect(self.ui.tab_b_ms2mzml_pb, QtCore.SIGNAL("clicked()"), self.b_load_mzml)
        QtCore.QObject.connect(self.ui.tab_b_saveimgfolder_pb, QtCore.SIGNAL("clicked()"), self.b_save_img2folder)
        QtCore.QObject.connect(self.ui.tab_b_sumxlsxpath_pb, QtCore.SIGNAL("clicked()"), self.b_save_output)
        QtCore.QObject.connect(self.ui.tab_b_cfgpath_pb, QtCore.SIGNAL("clicked()"), self.b_save_cfg)
        QtCore.QObject.connect(self.ui.tab_b_runhunter_pb, QtCore.SIGNAL("clicked()"), self.b_run_hunter)
        QtCore.QObject.connect(self.ui.tab_b_runcfg_pb, QtCore.SIGNAL("clicked()"), self.b_create_cfg)

        # slots for tab c
        QtCore.QObject.connect(self.ui.tab_c_mod_lst_pb, QtCore.SIGNAL("clicked()"), self.set_general_mod)
        QtCore.QObject.connect(self.ui.tab_c_fa_lst_pb, QtCore.SIGNAL("clicked()"), self.set_fa_white_lst)
        QtCore.QObject.connect(self.ui.tab_c_prostane_mod_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_mod)
        QtCore.QObject.connect(self.ui.tab_c_prostane_abbr_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_abbr)
        QtCore.QObject.connect(self.ui.tab_c_frag_pattern_pb, QtCore.SIGNAL("clicked()"), self.set_frag_pattern)
        QtCore.QObject.connect(self.ui.tab_c_hgcfg_pb, QtCore.SIGNAL("clicked()"), self.set_hg_specifc_pattern)
        QtCore.QObject.connect(self.ui.tab_c_scorecfg_pb, QtCore.SIGNAL("clicked()"), self.set_wfrag)
        QtCore.QObject.connect(self.ui.set_default_pb, QtCore.SIGNAL("clicked()"), self.set_default_cfg)
        QtCore.QObject.connect(self.ui.tab_c_calcsnr_pb, QtCore.SIGNAL("clicked()"), self.c_calc_snr_amp)

        # slots for tab d
        self.ui.tab_d_mutlimode_cmb.currentIndexChanged['QString'].connect(self.d_set_multi_mode)
        QtCore.QObject.connect(self.ui.tab_d_addcfg_pb, QtCore.SIGNAL("clicked()"), self.d_load_batchcfg)
        QtCore.QObject.connect(self.ui.tab_d_addcfgfolder_pb, QtCore.SIGNAL("clicked()"), self.d_load_batchcfgfolder)
        QtCore.QObject.connect(self.ui.tab_d_clearall_pb, QtCore.SIGNAL("clicked()"), self.ui.tab_d_infiles_pte.clear)
        QtCore.QObject.connect(self.ui.tab_d_runbatch_pb, QtCore.SIGNAL("clicked()"), self.d_run_batchmode)

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

    # def set_prostane_t(self):
    #     self.ui.label_12.show()
    #     self.ui.prostane_ox_yes_rb.show()
    #     self.ui.prostane_ox_no_rb.show()
    #
    # def set_prostane_f(self):
    #     self.ui.label_12.hide()
    #     self.ui.prostane_ox_yes_rb.hide()
    #     self.ui.prostane_ox_no_rb.hide()

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
        # if self.ui.prostane_ox_yes_rb.isChecked():
        #     ox_prostane_mode = 1
        # elif self.ui.prostane_ox_no_rb.isChecked():
        #     ox_prostane_mode = 0

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

    def b_save_cfg(self):
        b_save_cfg_path = QtGui.QFileDialog.getSaveFileName(caption='Save file', filter='.txt')
        self.ui.tab_b_cfgpath_le.clear()
        b_save_cfg_str = os.path.abspath(b_save_cfg_path[0])
        self.ui.tab_b_cfgpath_le.setText(b_save_cfg_str)

    def set_ms1_max(self):
        if self.ui.tab_b_ms1max_chb.isChecked():
            self.ui.tab_b_ms1max_lb.show()
            self.ui.tab_b_ms1max_spb.show()
        else:
            self.ui.tab_b_ms1max_lb.hide()
            self.ui.tab_b_ms1max_spb.hide()

    def b_get_params(self):
        # QtGui.QMessageBox.information(self, 'Information',
        #                               'LPPtiger need time to run, it might take up to few hours.\n'
        #                               'The Please click "OK" to start the analysis.')

        # self.ui.tab_b_loadlpppath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_lipids\Lv1_max3O_max1Keto_prostane\PC_Lv1_max3O_max1Keto_prostane.xlsx')
        # self.ui.tab_b_loadfapath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_lipids\Lv1_max3O_max1Keto_prostane\PC_Lv1_max3O_max1Keto_prostane_FA_SUM.xlsx')
        # self.ui.tab_b_ms2mzml_le.setText(r'D:\project_mzML\CM_DDA_neg_mzML\070120_CM_neg_70min_SIN_I.mzML')
        # self.ui.tab_b_saveimgfolder_le.setText(r'D:\Project_lpptiger\output_sdf\CM_LPPs\PC')
        # self.ui.tab_b_sumxlsxpath_le.setText(r'D:\Project_lpptiger\output_sdf\CM_LPPs\PC\CM_PC_Lv12.xlsx')
        # self.ui.tab_b_rtstart_dspb.setValue(24.0)
        # self.ui.tab_b_rtend_dspb.setValue(25.0)
        # self.ui.tab_b_msppm_spb.setValue(20)
        # self.ui.tab_b_ms2ppm_spb.setValue(50)
        # self.ui.tab_b_hgppm_spb.setValue(50)
        # self.ui.tab_b_dda_spb.setValue(12)
        # self.ui.tab_b_msthreshold_spb.setValue(2000)
        # self.ui.tab_b_ms2threshold_spb.setValue(1)
        # self.ui.tab_b_mzstart_dspb.setValue(850.0)
        # self.ui.tab_b_mzend_dspb.setValue(880.0)
        # self.ui.tab_b_ms2infoth_dspb.setValue(0)
        # self.ui.tab_b_ms2hginfoth_dspb.setValue(0)
        # self.ui.tab_c_snratio_spb.setValue(5)
        # self.ui.tab_c_ram_spb.setValue(16)

        if self.ui.tab_b_vendor_cmb.currentText() == 'Waters':
            usr_vendor = 'waters'
        elif self.ui.tab_b_vendor_cmb.currentText() == 'Thermo':
            usr_vendor = 'thermo'
        else:
            usr_vendor = ''
        if self.ui.mode_lcms_rb.isChecked():
            usr_exp_mode = 'LC-MS'
        elif self.ui.mode_static_rb.isChecked():
            usr_exp_mode = 'Shotgun'
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
        mzml_path_str = str(self.ui.tab_b_ms2mzml_le.text())
        img_output_folder_str = str(self.ui.tab_b_saveimgfolder_le.text())
        xlsx_output_path_str = str(self.ui.tab_b_sumxlsxpath_le.text())

        # tab_b_params
        rt_start = self.ui.tab_b_rtstart_dspb.value()
        rt_end = self.ui.tab_b_rtend_dspb.value()
        mz_start = self.ui.tab_b_mzstart_dspb.value()
        mz_end = self.ui.tab_b_mzend_dspb.value()
        dda_top = self.ui.tab_b_dda_spb.value()
        ms_th = self.ui.tab_b_msthreshold_spb.value()
        if self.ui.tab_b_ms1max_chb.isChecked():
            ms_max = self.ui.tab_b_ms1max_spb.value()
            print('Set ms1_max', ms_max)
        else:
            ms_max = 0
            print('No ms1_max', ms_max)
        ms2_th = self.ui.tab_b_ms2threshold_spb.value()
        ms_ppm = self.ui.tab_b_msppm_spb.value()
        ms2_ppm = self.ui.tab_b_ms2ppm_spb.value()
        pr_window = self.ui.tab_b_prwindow_spb.value()
        hg_th = self.ui.tab_b_hgthreshold_spb.value()
        hg_ppm = self.ui.tab_b_hgppm_spb.value()
        ms2_info_threshold = self.ui.tab_b_ms2infoth_dspb.value() * 0.01
        hgms2_info_threshold = self.ui.tab_b_ms2hginfoth_dspb.value() * 0.01
        overall_score_filter = self.ui.tab_b_score_spb.value()
        # from settings tab
        lipid_specific_cfg = self.ui.tab_c_hgcfg_le.text()
        score_cfg = self.ui.tab_c_scorecfg_le.text()
        sn_ratio = self.ui.tab_c_snratio_spb.value()
        isotope_score_filter = self.ui.tab_c_iso_score_spb.value()
        rank_score_filter = self.ui.tab_c_rank_score_spb.value()
        msp_score_filter = self.ui.tab_c_msp_score_spb.value()
        fp_score_filter = self.ui.tab_c_fp_score_spb.value()
        snr_score_filter = self.ui.tab_c_snr_score_spb.value()
        core_num = self.ui.tab_c_cores_spb.value()
        max_ram = self.ui.tab_c_ram_spb.value()
        img_typ = self.ui.tab_c_imagetype_cmb.currentText()[1:]
        img_dpi = self.ui.tab_c_dpi_spb.value()

        usr_parallelization_mode = self.ui.tab_c_parallization_cmb.currentIndex()
        if usr_parallelization_mode == 0:
            print('parallelization_mode: ', self.ui.tab_c_parallization_cmb.currentText())
            parallelization_mode = 'cpu'
        elif usr_parallelization_mode == 1:
            print('parallelization_mode: ', self.ui.tab_c_parallization_cmb.currentText())
            parallelization_mode = 'gpu'
        else:
            print('parallization_mode: ', self.ui.tab_c_parallization_cmb.currentText())
            parallelization_mode = 'cpu'

        usr_isotope_score_mode = self.ui.tab_c_isotopescoremode_cmb.currentIndex()
        if usr_isotope_score_mode == 0:
            print(self.ui.tab_c_isotopescoremode_cmb.currentText())
            fast_isotope = False
        else:
            print(self.ui.tab_c_isotopescoremode_cmb.currentText())
            fast_isotope = True

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

        start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        hunter_param_dct = {'hunter_folder': self.theolpp_cwd, 'hunter_start_time': start_time_str,
                            'vendor': usr_vendor, 'experiment_mode': usr_exp_mode,
                            'lipid_type': _pl_class, 'charge_mode': _pl_charge,
                            'lpp_sum_info_path_str': lpp_sum_info_path_str, 'fa_sum_path_str': fa_sum_path_str,
                            'mzml_path_str': mzml_path_str, 'img_output_folder_str': img_output_folder_str,
                            'xlsx_output_path_str': xlsx_output_path_str,
                            'rt_start': rt_start, 'rt_end': rt_end, 'mz_start': mz_start, 'mz_end': mz_end,
                            'ms_th': ms_th, 'ms_ppm': ms_ppm, 'pr_window': pr_window, 'dda_top': dda_top,
                            'ms2_th': ms2_th, 'ms2_ppm': ms2_ppm, 'ms2_infopeak_threshold': ms2_info_threshold,
                            'hg_th': hg_th, 'hg_ppm': hg_ppm, 'ms2_hginfopeak_threshold': hgms2_info_threshold,
                            'score_filter': overall_score_filter,
                            'lipid_specific_cfg': lipid_specific_cfg, 'score_cfg': score_cfg, 'sn_ratio': sn_ratio,
                            'isotope_score_filter': isotope_score_filter, 'rank_score_filter': rank_score_filter,
                            'msp_score_filter': msp_score_filter, 'fp_score_filter': fp_score_filter,
                            'snr_score_filter': snr_score_filter,
                            'parallization_mode': parallelization_mode, 'core_number': core_num, 'max_ram': max_ram,
                            'img_type': img_typ, 'img_dpi': img_dpi, 'fast_isotope': fast_isotope, 'ms_max': ms_max
                            }

        return hunter_param_dct

    def b_run_hunter(self):

        hunter_param_dct = self.b_get_params()

        print(hunter_param_dct)

        start_time = hunter_param_dct['hunter_start_time']

        param_log_output_path_str = (str(self.ui.tab_b_saveimgfolder_le.text()) +
                                     r'/LPPtiger_Params-Log_%s.txt' % start_time
                                     )
        if hunter_param_dct['vendor'] != '':
            config = configparser.ConfigParser()
            with open(param_log_output_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in hunter_param_dct.keys():
                    config.set('parameters', param, str(hunter_param_dct[param]))
                config.write(usr_param_cfg)

            tot_run_time = huntlipids(hunter_param_dct)

            if isinstance(tot_run_time, str):
                self.ui.tab_b_statusrun_pte.insertPlainText(tot_run_time)
                self.ui.tab_b_statusrun_pte.insertPlainText('>>> >>> >>> FINISHED <<< <<< <<<')

            else:
                self.ui.tab_b_statusrun_pte.insertPlainText('!! Failed read input files, '
                                                            'please check vendor settings!!')

        else:
            self.ui.tab_b_statusrun_pte.insertPlainText('Please check your vendor settings!')

    def b_create_cfg(self):

        hunter_param_dct = self.b_get_params()

        print(hunter_param_dct)

        param_log_output_path_str = str(self.ui.tab_b_cfgpath_le.text())
        if hunter_param_dct['vendor'] != '':
            config = configparser.ConfigParser()
            with open(param_log_output_path_str, 'w') as usr_param_cfg:
                config.add_section('parameters')
                for param in hunter_param_dct.keys():
                    config.set('parameters', param, str(hunter_param_dct[param]))
                config.write(usr_param_cfg)
                self.ui.tab_b_statuscfg_pte.insertPlainText('>>> >>> >>> SAVED <<< <<< <<<')
        else:
            self.ui.tab_b_statuscfg_pte.insertPlainText('Please check your vendor settings!')

    def c_calc_snr_amp(self):
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

    def d_set_multi_mode(self):

        multi_mode_idx = self.ui.tab_d_mutlimode_cmb.currentIndex()

        if multi_mode_idx == 1:
            print('Multi processing mode')
            self.ui.tab_d_maxbatch_lb.show()
            self.ui.tab_d_maxbatch_spb.show()
            self.ui.tab_d_maxsubcore_lb.show()
            self.ui.tab_d_maxsubcore_spb.show()
            self.ui.tab_d_maxsubram_lb.show()
            self.ui.tab_d_maxsubram_spb.show()
        elif multi_mode_idx == 0:
            print('Single processing mode')
            self.ui.tab_d_maxbatch_lb.hide()
            self.ui.tab_d_maxbatch_spb.hide()
            self.ui.tab_d_maxsubcore_lb.hide()
            self.ui.tab_d_maxsubcore_spb.hide()
            self.ui.tab_d_maxsubram_lb.hide()
            self.ui.tab_d_maxsubram_spb.hide()

    @staticmethod
    def d_get_same_files(folder, filetype_lst):
        """
        find all files with same type in specified folder
        :param str folder: absolute file path
        :param list filetype_lst: e.g. ['*.mzml', '*.mzML']
        :return: a list of absolute file path
        :rtype: list
        """
        if folder is not u'':
            os.chdir(folder)
            _pre_found_lst = []
            for _filetype in filetype_lst:
                _tmp_found_lst = glob.glob(_filetype)
                # merge list
                _pre_found_lst += [f for f in _tmp_found_lst if f not in _pre_found_lst]
            filename_lst = _pre_found_lst
            abs_path_lst = list(os.path.abspath(ff) for ff in _pre_found_lst)
        else:
            filename_lst = []
            abs_path_lst = []

        return filename_lst, abs_path_lst

    def d_load_batchcfg(self):
        # check existed files
        _loaded_files = str(self.ui.tab_d_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfg_dialog = QtGui.QFileDialog(self)
        b_load_cfg_dialog.setNameFilters([u'LPPtiger batch mode files (*.txt)'])
        b_load_cfg_dialog.selectNameFilter(u'LPPtiger batch mode files (*.txt)')
        if b_load_cfg_dialog.exec_():
            b_load_cfg_str = b_load_cfg_dialog.selectedFiles()[0]
            b_load_cfg_str = os.path.abspath(b_load_cfg_str)
            if b_load_cfg_str not in _loaded_lst:
                self.ui.tab_d_infiles_pte.insertPlainText(b_load_cfg_str)  # take unicode only
                self.ui.tab_d_infiles_pte.insertPlainText(u'\n')
            else:
                _msgBox = QtGui.QMessageBox()
                _msgBox.setText(u'Batch config file has been chosen already.')
                _msgBox.exec_()

    def d_load_batchcfgfolder(self):
        # check existed files
        _loaded_files = str(self.ui.tab_d_infiles_pte.toPlainText())
        _loaded_lst = _loaded_files.split('\n')

        b_load_cfgfolder_str = QtGui.QFileDialog.getExistingDirectory()
        _cfg_name_lst, _cfg_path_lst = self.d_get_same_files(b_load_cfgfolder_str, filetype_lst=['*.txt', '*.txt'])
        _duplicated_str = ''
        for _cfg in _cfg_path_lst:
            if _cfg not in _loaded_lst:
                self.ui.tab_d_infiles_pte.insertPlainText(_cfg)
                self.ui.tab_d_infiles_pte.insertPlainText('\n')
            else:
                _duplicated_str = _duplicated_str + _cfg + '\n'
        if len(_duplicated_str) > 0:
            _msgBox = QtGui.QMessageBox()
            _msgBox.setText(_duplicated_str + u'Already chosen. \n Skipped')
            _msgBox.exec_()

    @staticmethod
    def d_read_batch_cfg(batch_cfg):
        config = configparser.ConfigParser()
        config.read(batch_cfg)
        batch_cfg_dct = {}
        if config.has_section('parameters'):
            user_cfg = 'parameters'
            options = config.options(user_cfg)
            for param in options:
                batch_cfg_dct[param] = config.get(user_cfg, param)
        batch_cfg_key_lst = batch_cfg_dct.keys()
        i_type_key_lst = ['ms_th', 'ms2_th', 'hg_th', 'ms_ppm', 'ms2_ppm', 'hg_ppm', 'dda_top', 'sn_ratio',
                          'core_number', 'max_ram', 'img_dpi', 'ms_max']
        f_type_key_lst = ['rt_start', 'rt_end', 'mz_start', 'mz_end', 'pr_window',
                          'ms2_infopeak_threshold', 'ms2_hginfopeak_threshold',
                          'score_filter', 'isotope_score_filter', 'rank_score_filter',
                          'msp_score_filter', 'fp_score_filter', 'snr_score_filter']

        if len(batch_cfg_key_lst) > 0:
            for cfg_key in batch_cfg_key_lst:
                if cfg_key in i_type_key_lst:
                    try:
                        batch_cfg_dct[cfg_key] = int(batch_cfg_dct[cfg_key])
                    except ValueError:
                        batch_cfg_dct[cfg_key] = int(float(batch_cfg_dct[cfg_key]))
                elif cfg_key in f_type_key_lst:
                    batch_cfg_dct[cfg_key] = float(batch_cfg_dct[cfg_key])

        return batch_cfg_dct

    def d_run_batchmode(self):

        self.ui.tab_d_statusrun_pte.clear()

        loaded_cfg_files = str(self.ui.tab_d_infiles_pte.toPlainText())
        pre_loaded_cfg_lst = loaded_cfg_files.split('\n')

        loaded_cfg_lst = []
        for f in pre_loaded_cfg_lst:
            if len(f) > 4:
                loaded_cfg_lst.append(f)

        tot_num = len(loaded_cfg_lst)
        run_counter = 1

        multi_mode_idx = self.ui.tab_d_mutlimode_cmb.currentIndex()

        os.chdir(self.theolpp_cwd)

        if multi_mode_idx == 1:  # multi mode

            max_process = self.ui.tab_d_maxbatch_spb.value()
            sub_max_core = self.ui.tab_d_maxsubcore_spb.value()
            sub_max_ram = self.ui.tab_d_maxsubram_spb.value()

            cfg_dct_lst = []
            for cfg_file in loaded_cfg_lst:
                hunter_param_dct = self.d_read_batch_cfg(cfg_file)
                if 'vendor' in hunter_param_dct.keys():
                    hunter_param_dct['batch_cfg_file'] = cfg_file
                    hunter_param_dct['core_number'] = sub_max_core
                    hunter_param_dct['max_ram'] = sub_max_ram
                    cfg_dct_lst.append(hunter_param_dct)
                else:
                    hunter_param_dct['batch_cfg_file'] = ''

            if len(cfg_dct_lst) > max_process:
                sub_part_lst = map(None, *(iter(cfg_dct_lst),) * max_process)
            else:
                sub_part_lst = [cfg_dct_lst]

            tot_part = len(sub_part_lst)
            part_num = 1
            for sub_cfg_lst in sub_part_lst:
                sub_cfg_lst = filter(lambda x: x is not None, sub_cfg_lst)
                parallel_pool = multiprocessing.pool.ThreadPool(max_process)
                hunter_results_lst = []
                core_worker_count = 1
                for _cfg_dct in sub_cfg_lst:
                    time.sleep(1)
                    self.ui.tab_d_statusrun_pte.insertPlainText('Start Batch %i / %i file %i / %i ...\n'
                                                                % (part_num, tot_part, core_worker_count, max_process))
                    self.ui.tab_d_statusrun_pte.insertPlainText('>>> processing...\n')

                    start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                    _cfg_dct['hunter_start_time'] = start_time_str
                    os.chdir(_cfg_dct['hunter_folder'])
                    tot_run_time = parallel_pool.apply_async(huntlipids, args=(_cfg_dct,))

                    core_worker_count += 1
                    hunter_results_lst.append(tot_run_time)

                parallel_pool.close()
                parallel_pool.join()

                for hunter_time in hunter_results_lst:

                    run_time = str(hunter_time.get())

                    if isinstance(run_time, str):
                        self.ui.tab_d_statusrun_pte.appendPlainText('>>> %s' % run_time)
                        self.ui.tab_d_statusrun_pte.appendPlainText('FINISHED with file %i / %i\n' %
                                                                    (run_counter, tot_num))
                        run_counter += 1
                    else:
                        self.ui.tab_d_statusrun_pte.insertPlainText(
                            '!! Failed to process batch mode configure file:\n Please check settings!!')

                part_num += 1

        else:  # single mode
            for _cfg in loaded_cfg_lst:

                self.ui.tab_d_statusrun_pte.insertPlainText('Start processing...\n%s\n' % _cfg)
                hunter_param_dct = self.d_read_batch_cfg(_cfg)
                if 'vendor' in hunter_param_dct.keys():
                    start_time_str = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
                    hunter_param_dct['hunter_start_time'] = start_time_str
                    tot_run_time = huntlipids(hunter_param_dct)

                    if isinstance(tot_run_time, str):
                        self.ui.tab_d_statusrun_pte.appendPlainText(tot_run_time)
                        self.ui.tab_d_statusrun_pte.appendPlainText('FINISHED with file %i / %i\n' %
                                                                    (run_counter, tot_num))
                        run_counter += 1
                    else:
                        self.ui.tab_d_statusrun_pte.insertPlainText(
                            '!! Failed read batch mode configure files:\n %s \n Please check settings!!' % _cfg)
                else:
                    self.ui.tab_d_statusrun_pte.insertPlainText(
                        '!! Failed read batch mode configure files:\n %s \n Please check settings!!' % _cfg)


if __name__ == '__main__':
    multiprocessing.freeze_support()
    import sys
    app = QtGui.QApplication(sys.argv)
    window = LPPtigerMain()
    window.show()
    sys.exit(app.exec_())
