# -*- coding: utf-8 -*-
# Copyright 2015 Zhixu Ni, AG Bioanalytik, BBZ, University of Leipzig.
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release of oxLPPdb.
# For more info please contact: zhixu.ni@uni-leipzig.de

import re


class DBoxTheo:
    def __init__(self):
        print("Start to oxidize-->")
        self.theo_ox_dct = {}

    @staticmethod
    def hydroxyl():
        hydroxyl_lst = ['C(O)C']

    @staticmethod
    def release():
        print("  locker.release() called.|no object needed|")


def bulk_oxidizer(cls):
    """

    :param cls:
    :return:
    """

    def _oxidizer(func):
        def __oxidizer():
            print("before %s called [%s]." % (func.__name__, cls))
            cls.acquire()
            try:
                return func()
            finally:
                cls.release()

        return __oxidizer

    return _oxidizer


@oxidizer(DBoxTheo)
def hydroxyl():
    print(" myfunc() called.")


hydroxyl()
hydroxyl()
