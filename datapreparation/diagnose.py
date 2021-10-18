#! /usr/bin/env python3


""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2013-03-26 version

    bonsai/src_v02/diagnose.py
    processing the diagnosis data

    Notice: This file is not imported
    using the name dia, since dia is
    often used for a dataframe with
    diagnosis data content

------------------------------------"""

import pandas as pd
import numpy as np
import copy

import bonsai_io as bio
import common
import lexicon
import global_settings as gs
import merge





""" ----------------------------------------

    generate version 1 dia 

---------------------------------------"""



def special(ICD10):
    if common.isnull(ICD10): # not needed after filled in
        return 0
    if 'G' in str(ICD10):
        return 1    
    if 'H' in str(ICD10):
        return 1    
    if 'K' in str(ICD10):
        return 1 
    if 'L' in str(ICD10):
        return 1
    if 'O' in str(ICD10):
        return 1
    if 'Q' in str(ICD10):
        return 1
    if 'R' in str(ICD10):
        return 1
    if 'T' in str(ICD10):
        return 1
    if 'Z' in str(ICD10):
        return 1
    if str(ICD10) == 'D469':
        return 0
    if str(ICD10) == 'D761':
        return 0
    if str(ICD10) == 'D459':
        return 0
    if 'D' in str(ICD10):
        return 1
    return 0


def generate_dia():
    """
    constructing the file stored as generated_dia,
    see places.py
    """

    xcols  =    ['ICD7', 'ICD9', 'text_C24_']

    dia = bio.read_original_dia()
    dia = dia.sort_values(by = xcols)
    print('orig shape:', dia.shape)

    # (1) select the unique rows
    dia = dia.drop_duplicates()
    print('non duplicate shape:', dia.shape)

    # (2) select first diagnoses
    dia = dia[dia['DiagnosNr_Diagnos'] == '1']
    dia = dia.drop(['DiagnosNr_Diagnos'], axis=1)
    print('first dia shape:',dia.shape)

    # (3) fill in codes
    dia = lexicon.fill_in_by_compiled_lex(dia)
    dia = dia.drop(xcols, axis=1)
    print('filled dia shape:',dia.shape)
    
    # (4) remove the special cases and the not needed columns
    dia['special'] = dia['ICD10'].apply(special)
    dia = dia[dia['special'] == 0]
    dia = dia.drop(['special'], axis=1)
    print('no special shape:',dia.shape)

    # (5) take care of numbers
    if 'Diagnos_lder' in gs.places.diagnose_selection:
        dia['Diagnos_lder'] = dia['Diagnos_lder'].apply(common.str2number)


    return dia


def rm_dia_cols(dia):
    cols = gs.places.diagnose_ohe
    dia = common.rmcols(dia, cols)
    return dia



""" ----------------------------------------

    add dia groups

---------------------------------------"""



def look_up_group_by_codes(groups, ICD10, SNOMED):
    g = gs.names.input_data_group
    row = groups[(groups['ICD10'] == ICD10)  & (groups['SNOMED'] == SNOMED)]
    if not row.empty:
        return row[g].values[0]
    return str(0)

def look_up_group(df, groups):
    g = look_up_group_by_codes(groups, df['ICD10'], df['SNOMED'])
    return g

def add_group(dia, groups):
    g   = gs.names.dia_group
    dia[g] = dia.apply(lambda d: look_up_group(d, groups), axis = 1)
    return dia

def rm_group_col(dia):
    g   = gs.names.dia_group
    dia = common.rmcol(dia, g)
    return dia


""" ----------------------------------------

    to one data frame df1 add columns 
    from another data frame df2

    Note: in add_cols
       df1[yc] = 0    when there is no entry x in df2
       df1[yc] = NaN  when df2[x] = NaN 

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



def look_up_or_zero(entry, df, entry_col, value_col, verb = False):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        val = dfe[value_col].values[0]
        if verb: 
            print('val =', val, type(val))
        if isinstance(val, str):
            return val
    return str(0)

def add_cols_or_zero(df1, df2, xc, ycs, verb = False):
    df1_copy = copy.copy(df1) 
    for yc in ycs:
        df1_copy[yc] = df1[xc].apply(lambda x:look_up_or_zero(x, df2, xc, yc, verb))
    return df1_copy

"""
def add_cols_or_zero(df1, df2, xc, ycs, verb = False):
    for yc in ycs:
        df1[yc] = df1[xc].apply(lambda x:look_up_or_zero(x, df2, xc, yc, verb))
    return df1
"""

def look_up_aho(df, x, col, zero = False):
    return df[col][x] if x in df.index and (not zero or isinstance(df[col][x], str)) else str(0)

def add_cols_aho(df1, df2, xc, ycs, zero = False):
    if zero:
        for yc in ycs:
            df1[yc] = df1[xc].apply(lambda x: df2[yc][x] if x in df2.index and isinstance(df2[yc][x], str) else str(0))
    else:
        for yc in ycs:
            df1[yc] = df1[xc].apply(lambda x: df2[yc][x] if x in df2.index else str(0))
    return df1
    
""" ----------------------------------------

    add person data

---------------------------------------"""


def add_pers(dia):
    cols = gs.places.person_cols  
    copy_cols = copy.copy(cols)    
    copy_cols.remove('LopNr')
    pers = bio.readperson()
    dia = add_cols(dia, pers, 'LopNr', copy_cols)
    return dia

def add_pers_ohe(dia):
    cols = gs.places.person_ohe
    df = ones_x(dia, cols)
    return df

def rm_pers_cols(dia):
    cols = gs.places.person_ohe
    dia = common.rmcols(dia, cols)
    return dia



""" ----------------------------------------

    add incare data (sluten vaard)

---------------------------------------"""

def add_incare(dia, nr = 0):

    xcol = 'LopNr'
    ycol = gs.places.incare_ohe[0]          # only one column is used so far

    inc = bio.readincare()
    inc = inc.sort_values([xcol]).reset_index(drop=True)
    if nr > 0:
        inc = inc[0:nr]                     # an initial part of incare
    L = list(dia[xcol].values)
    inc = inc[ inc[xcol].isin(L) ]          # part of inc with LopNr in dia
    inc = name_compression(inc, xcol, ycol) # first letter set for each LopNr
    dia = add_cols(dia, inc, xcol, [ycol])  # add compressed inc cols to dia    
    return dia

#   the following functions are just for test since the incare compression
#   lists are merged with the nicare and causes lists before unfolding ohe


def add_incare_ohe(dia):
    ycol = gs.places.incare_ohe[0]          # only one column is used so far
    dia  = one_general(dia, ycol)           # mk first letter one hot
    return dia

def rm_incare_cols(dia):
    cols = gs.places.incare_ohe
    dia = common.rmcols(dia, cols)
    return dia



""" ----------------------------------------

    add nicare data (oppen vaard)

---------------------------------------"""


def add_nicare(dia, nr = 0):

    xcol = 'LopNr'
    ycol = gs.places.nicare_ohe[0]          # only one column is used so far

    nic = bio.readnicare()
    nic = nic.sort_values([xcol]).reset_index(drop=True)
    if nr > 0:
        nic = nic[0:nr]                     # an initial part of incare
    L = list(dia[xcol].values)
    nic = nic[ nic[xcol].isin(L) ]          # part of nic with LopNr in dia
    nic = name_compression(nic, xcol, ycol) # first letter set for each LopNr
    dia = add_cols(dia, nic, xcol, [ycol])  # add compressed nic cols to dia    
    return dia

#   the following functions are just for test since the nicare compression
#   lists are merged with the nicare and causes lists before unfolding ohe

def add_nicare_ohe(dia):
    ycol = gs.places.nicare_ohe[0]          # only one column is used so far
    dia  = one_general(dia, ycol)           # mk first letter one hot
    return dia

def rm_nicare_cols(dia):
    cols = gs.places.nicare_ohe
    dia = common.rmcols(dia, cols)
    return dia


""" ----------------------------------------

    add drug data

---------------------------------------"""


def add_drug(dia):
    cols = gs.places.drug_selection.copy()
    cols.remove('LopNr')
    drug = bio.readdrug()
    dia = add_cols(dia, drug, 'LopNr', cols)
    return dia


def add_drug_ohe(dia):
    cols = gs.places.drug_ohe
    df = ones_x(dia, cols)
    return df


def rm_drug_cols(dia):
    cols = gs.places.drug_ohe
    dia = common.rmcols(dia, cols)
    return dia



""" ----------------------------------------

    add one hot encodings for names column
    with unique names (LopNr) and a single
    code in each row in the codes column

---------------------------------------"""

"""

def equal_str(a, b):
    if not (isinstance(a, str) and isinstance(b, str)):
        return 0
    return int(a == b)

"""


def one(df, c):
    for x in  df[c].dropna().drop_duplicates():
        df[x] = (df[c] == x).astype(int)
    return df

"""
def one_x(df, c):
    for x in  df[c].drop_duplicates():
        if isinstance(x, str):
            df[c + '_' + x] = df[c].apply(lambda z: equal_str(z, x)) 
    return df
"""

def one_x(df, c):
    for x in  df[c].dropna().drop_duplicates():
        df[c + '_' + x] = (df[c] == x).astype(int)
    return df

    

def ones_x(df, cs):
    for c in cs:
        df = one_x(df, c)
    return df


def to_int(x):
    if isinstance(x, str):
        return int(x)
    if isinstance(x, int):
        return x
    return -1

def nr_sort(xs):
    """
    sort a list of numbers on str type
    """
    if not common.isarray(xs):
        return []
    ixs = list(map(to_int, xs))
    ixs.sort()
    xs = list(map(str, ixs))
    return xs

def one_sorted(df, c):
    xs = list(df[c].drop_duplicates())
    xs = nr_sort(xs)
    for x in xs:
        df[c + '_' + x] = (df[c] == x).astype(int)
    return df
    
def add_one_hot_groups(dia):
    grp = gs.names.dia_group
    dia = one_sorted(dia, grp)
    return dia



""" ----------------------------------------

    add one hot encodings for names column with 
    non-unique names and possibly several space
    separated codes in each row in the codes 
    column

---------------------------------------"""

def head(xs):
    ys = []
    for x in xs:
        if isinstance(x, str):
            ys += [x[0]]
    return ys

def split_and_head(ys):
    cs = []
    for y in ys.values:
        if not common.isnull(y):
            cs = np.append(cs, y.split())
    return np.unique(head(cs))


def name_compression(df, xcol, ycol):
    data = []
    xs =  df[xcol].drop_duplicates()
    for x in xs:
        dx = df[ df[xcol] == x ]        
        ys = dx[ycol].drop_duplicates()
        ys = split_and_head(ys)
        data = data + [ [x] + [ys] ]
    ds = pd.DataFrame(data)
    ds.columns = [xcol, ycol]    
    return ds



def unique_values(ys):
    ys1 = [] 
    for y in ys.values:
        ys1 = np.append(ys1, np.array(y))
    return np.unique(ys1)

def in_col_list(d, col, y):
    if isinstance(d[col], str):
        return 0
    if common.isarray(d[col]):
        return int(y in list(d[col]))
    print('ERROR ? in_col_list: d[col] is neither a str or an array')
    print('d[col] =', d[col], 'type =', type(d[col]),'\ty=', y)
    return 0

def one_general(df, ycol):
    # f  = lambda d: ( int(y in d[ycol]))
    ys = unique_values(df[ycol])
    for y in ys:
        df[ycol + '_' +y] = df.apply(lambda d:in_col_list(d, ycol, y), axis = 1)
    return df
