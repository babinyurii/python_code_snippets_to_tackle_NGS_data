#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 15:58:56 2018
@author: yuriy
"""
from math import log2

def log2_pseudo_count(arr):
    """returns pd.Series obj of log2 
    transformed values where 0.0 are substituted 
    with pseudocount (minimum value in the Series/10)
    
    arr is a pandas Series obj
    """
    
    min_value = arr[arr > 0].min()
    pseudo_count = min_value / 10
    arr_log2 = []
    for num in arr.values:
        if num == 0:
            arr_log2.append(log2(pseudo_count))
        else:
            arr_log2.append(log2(num))
    return pd.Series(arr_log2)



