#! /usr/bin/env python3


""" -------------------------------

    merge.py

    Copyright (C) 2018 RISE
    This code was produced by RISE

    First version created 2013-04-05
    based on code from other files

    bonsai/src_v02/merge.py

------------------------------------"""

import pandas as pd
import numpy as np

import bonsai_io as bio
import common

import global_settings as gs





""" ----------------------------------------

    to one data frame df1 add columns 
    from another data frame df2

---------------------------------------"""

def look_up_entry(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        return dfe[value_col].values[0]
    return str(0)

def add_cols(df1, df2, xc, ycs):
    for yc in ycs:
        df1[yc] = df1[xc].apply(lambda x:look_up_entry(x, df2, xc, yc))
    return df1


""" ----------------------------------------

    add person data

---------------------------------------"""


def add_pers(df):
    cols = gs.places.person_selection.copy()
    cols.remove('LopNr')
    pers = bio.readperson()
    df = add_cols(df, pers, 'LopNr', cols)
    # cols = gs.places.person_ohe
    return df

def add_pers_ohe(df):
    cols = gs.places.person_ohe
    df = ones_x(df, cols)
    return df

def rm_pers_cols(df):
    cols = gs.places.person_ohe
    df = common.rmcols(df, cols)
    return df

