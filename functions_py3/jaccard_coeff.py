# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:14:15 2019

@author: Meghana
"""


def jaccard_coeff(set1, set2):
    ls1 = len(set1)
    ls2 = len(set2)
    if ls1 == 0 and ls2 == 0:
        return 1
    inter = len(set1.intersection(set2))
    return float(inter) / (ls1 + ls2 - inter)
