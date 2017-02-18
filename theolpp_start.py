# -*- coding: utf-8 -*-
# Copyright 2016-2017 LPP team, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of LipidHunter.
# For more info please contact:
#     LPP team oxlpp@bbz.uni-leipzig.de
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de


import os
from PySide import QtGui
from LibTheoLPP.TheoLPP_Main import TheoLPP_Main
import sys

if __name__ == '__main__':
    usr_cwd = os.getcwd()
    gui = QtGui.QApplication(sys.argv)
    LipidHunter = TheoLPP_Main()
    LipidHunter.show()
    sys.exit(gui.exec_())
