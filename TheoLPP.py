# -*- coding: utf-8 -*-
# Copyright 2015-2016 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

try:  # python3
    import configparser
except NameError:  # python2
    import ConfigParser as configparser

import os

import xlrd
from PySide import QtCore, QtGui

from TheoLPP_UI import Ui_MainWindow
from TheoLPP_Core import theolpp


class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

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
        QtCore.QObject.connect(self.ui.mod_lst__pb, QtCore.SIGNAL("clicked()"), self.set_general_mod)
        QtCore.QObject.connect(self.ui.fa_lst_pb, QtCore.SIGNAL("clicked()"), self.set_fa_white_lst)
        QtCore.QObject.connect(self.ui.prostane_mod_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_mod)
        QtCore.QObject.connect(self.ui.prostane_abbr_lst_pb, QtCore.SIGNAL("clicked()"), self.set_prostane_abbr)
        QtCore.QObject.connect(self.ui.frag_pattern_pb, QtCore.SIGNAL("clicked()"), self.set_frag_pattern)
        QtCore.QObject.connect(self.ui.set_default_pb, QtCore.SIGNAL("clicked()"), self.set_default_cfg)

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
                self.ui.mod_lst_le.setText(config.get(user_cfg, 'general_mod_lst'))
            if 'fa_white_list' in options:
                self.ui.fa_lst_le.setText(config.get(user_cfg, 'fa_white_list'))
            if 'prostane_mod_lst' in options:
                self.ui.prostane_mod_lst_le.setText(config.get(user_cfg, 'prostane_mod_lst'))
            if 'prostane_abbr' in options:
                self.ui.prostane_abbr_lst_le.setText(config.get(user_cfg, 'prostane_abbr'))
            if 'frag_patterns' in options:
                self.ui.frag_pattern_le.setText(config.get(user_cfg, 'frag_patterns'))

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
        self.ui.mod_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.mod_lst_le.setText(a_load_table_str)

    def set_fa_white_lst(self):
        self.ui.fa_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.fa_lst_le.setText(a_load_table_str)

    def set_prostane_mod(self):
        self.ui.prostane_mod_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.prostane_mod_lst_le.setText(a_load_table_str)

    def set_prostane_abbr(self):
        self.ui.prostane_abbr_lst_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.prostane_abbr_lst_le.setText(a_load_table_str)

    def set_frag_pattern(self):
        self.ui.frag_pattern_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.csv *.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.csv *.xlsx)')
        if load_table_dialog.exec_():
            a_load_table_str = load_table_dialog.selectedFiles()[0]
            a_load_table_str = os.path.abspath(a_load_table_str)
            self.ui.frag_pattern_le.setText(a_load_table_str)

    def set_default_cfg(self):
        config = configparser.ConfigParser()
        with open('configure.ini', 'w') as default_cfg:
            config.add_section('settings')
            config.set('settings', 'general_mod_lst', self.ui.mod_lst_le.text())
            config.set('settings', 'fa_white_list', self.ui.fa_lst_le.text())
            config.set('settings', 'prostane_mod_lst', self.ui.prostane_mod_lst_le.text())
            config.set('settings', 'prostane_abbr', self.ui.prostane_abbr_lst_le.text())
            config.set('settings', 'frag_patterns', self.ui.frag_pattern_le.text())
            config.write(default_cfg)

    def load_lipid_list(self):
        self.ui.load_lipid_le.clear()
        load_table_dialog = QtGui.QFileDialog(self)
        load_table_dialog.setNameFilters([u'Tables (*.xlsx)'])
        load_table_dialog.selectNameFilter(u'Tables (*.xlsx)')
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
        save_sdf_path = QtGui.QFileDialog.getSaveFileName(caption=u'Save file', filter=u'.sdf')
        self.ui.save_sdf_le.clear()
        save_sdf_str = os.path.abspath(save_sdf_path[0])
        self.ui.save_sdf_le.setText(save_sdf_str)

    def save_msp(self):
        save_msp_path = QtGui.QFileDialog.getSaveFileName(caption=u'Save file', filter=u'.msp')
        self.ui.save_msp_le.clear()
        save_msp_str = os.path.abspath(save_msp_path[0])
        self.ui.save_msp_le.setText(save_msp_str)

    def run_theolpp(self):

        self.ui.run_status_te.append(u'start!')
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

        msp_mode = 0
        if self.ui.prostane_ox_yes_rb.isChecked():
            msp_mode = 1
        elif self.ui.prostane_ox_no_rb.isChecked():
            msp_mode = 0

        lipid_lst_path = self.ui.load_lipid_le.text()

        lipid_tab = self.ui.exceltab_cb.currentText()
        lipid_class = self.ui.pl_class_cb.currentText()
        lipid_class = lipid_class[-3:-1]  # get PL abbr

        sdf_path = self.ui.save_sdf_le.text()
        msp_path = self.ui.save_msp_le.text()
        mod_lst_path = self.ui.mod_lst_le.text()
        fa_lst_path = self.ui.fa_lst_le.text()
        prostane_mod_path = self.ui.prostane_mod_lst_le.text()
        prostane_abbr_path = self.ui.prostane_abbr_lst_le.text()
        frag_pattern_path = self.ui.frag_pattern_le.text()

        param_dct = {'lipid_class': lipid_class, 'ox_level': ox_level,
                     'oap_mode': oap_mode, 'ocp_mode': ocp_mode,
                     'lyso_oap_mode': lyso_oap_mode, 'lyso_ocp_mode': lyso_ocp_mode,
                     'ox_max': ox_max, 'keto_max': keto_max, 'ooh_max': ooh_max, 'epoxy_max': epoxy_max,
                     'lipid_lst_path': lipid_lst_path, 'lipid_tab': lipid_tab,
                     'prostane_mode': prostane_mode, 'ox_prostane_mode': ox_prostane_mode,
                     'sdf_path': sdf_path, 'msp_mode': msp_mode, 'msp_path': msp_path,
                     'mod_lst_path': mod_lst_path, 'fa_lst_path': fa_lst_path, 'prostane_mod_path': prostane_mod_path,
                     'prostane_abbr_path': prostane_abbr_path, 'frag_pattern_path': frag_pattern_path}

        # print(param_dct)
        info_1, info_2 = theolpp(param_dct)
        self.ui.run_status_te.append(info_1)
        self.ui.run_status_te.append(info_2)
        self.ui.run_status_te.append(u'Finished!')


if __name__ == '__main__':
    import sys

    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
