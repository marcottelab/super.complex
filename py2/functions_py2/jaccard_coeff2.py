# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:30:11 2020

@author: Meghana
"""

def jaccard_coeff2(list1, list2):
    inter = len(set(list1).intersection(set(list2)))
    union = len(list1) + len(list2) - inter
    return float(inter)/union