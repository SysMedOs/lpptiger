# -*- coding: utf-8 -*-
#
# Copyright 2016-2017 SysMedOs team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LPPtiger.
# For more info please contact:
#     SysMedOs team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#


import os
import sys

from PySide import QtGui

from LibLPPtiger.LPPtiger_Main import LPPtiger_Main

if __name__ == '__main__':
    usr_cwd = os.getcwd()
    gui = QtGui.QApplication(sys.argv)
    LPPtiger = LPPtiger_Main(cwd=usr_cwd)
    LPPtiger.show()
    sys.exit(gui.exec_())