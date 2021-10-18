#! /usr/bin/env python3


""" -------------------------------

    analyse.py

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2013-04-10 version

    bonsai/src_v02/analyze.py

    simple analysis of pandas dataframes data
    such as 

       1. find duplicated rows

       2. number of unique values in a column

       3. number of unique values in common 
          between two columns in two different
          files
       
       4. 
   
------------------------------------"""


import global_settings as gs
import numpy as np
import pandas as pd
import bonsai_io as bio
import common
import copy

def nr_of_unique_rows(df):
    d = df.drop_duplicates()
    return len(d)

def nr_of_unique_values_in_cols(df, cols):
    c = df.drop_duplicates(subset = cols)
    return len(c)


def nr_of_unique_values(df, col):
    c = df[col].dropna()
    c = c.drop_duplicates()
    return len(c)

"""
def nr_of_unique_numeric_values(df, col):

    c = df[col].dropna()
    c = c.drop_duplicates()
    c = c.str.isnumeric() 
    c = c[c].index.values
"""


def nr_of_nonnan_values(df, col):

    c = df[col].dropna()
    return len(c)
 
def nr_of_unique_digital_values(df, col):

    c = df[col].dropna()
    c = c.drop_duplicates()
    c = c.str.isdigit() 
    c = c[c].index.values
    # df = df.drop_duplicates(subset = col)
    # df = df[ df[col].dropna().str.isdigit() ]
    # df = df[ df[col].str.contains('\d', regex=True) ]
    return len(c)

def duplicated_rows(df):
    df['dup'] = df.duplicated()
    df = df[df['dup'] == True]
    return df

def print_duplicated_rows(df, nr):
    dup = duplicated_rows(df)
    print('Nr of rows in total', len(df))
    print('Nr of duplicated rows', len(dup))
    nr = min( nr,len(dup) )
    if nr > 0:
        print('the first', nr,' of them')
        print(dup[0:nr])
    return dup

def unique_number_values(df, col):
    df = df.drop_duplicates(subset = col)
    df = df[ df[col].str.contains('\d', regex=True) ]
    return df


def info(df, name = ''):
    print()
    if name != '':
        print()
        print('--------------------------------------------------')
        print()
        print('\tInfo on the file\n\t' + name)
        print()
        print('--------------------------------------------------')
        print()
        df_unique_nr = nr_of_unique_rows(df)
    print('          shape', df.shape)
    print('    unique rows', df_unique_nr)

    for c in df.columns:
        print()
        print('\tInfo on non-nan values of column', c)
        print()
        nonnan_nr = nr_of_nonnan_values(df, c)
        unique_nr = nr_of_unique_values(df, c)
        digital_nr =  nr_of_unique_digital_values(df, c)
        # numeric_nr = nr_of_unique_numeric_values(df, c)
        print('non-nan values', nonnan_nr)
        print(' unique values', unique_nr)
        print('digital values', digital_nr)
        # print('numeric values', unique_nr)
        
    print()
    # return unique_number_values(df, 'ICD10')

# df = df[ df[c].str.contains('\d', regex=True) ]



def readall():
    dia = bio.read_generated_dia()
    dgr = bio.read_diagroups()
    per = bio.readperson()
    ctr = bio.readcontrol()
    inc = bio.readincare()
    nic = bio.readnicare()
    dru = bio.readdrug()
    dcl = bio.readdrugclasses()
    tre = bio.readtreatment()
    sur = bio.readsurgery()
    cau = bio.readcause()

    data  = [
        dia, 
        dgr, 
        per,
        ctr, 
        inc, 
        nic, 
        dru, 
        dcl, 
        tre,
        sur,
        cau
]

    name = [
        'diagnos          ',
        'diagnosgrupp     ',
        'person           ',
        'kontrollgrupp    ',
        'sluten v_rd      ',
        '_ppen v_rd       ',
        'l_kemedel        ',
        'l_kemedelsgrupper',
        'behandling       ',
        'kirurgi          ',
        'orsak            ',
    ]

    return data, name


def info_on_all():

    data, name = readall()
        
    for i in range(0, len(name)):
        info(data[i], name[i])


def compare_lopnr(dfx, dfy, namex = 'data 1', namey = 'data 2'):

    xs = list(dfx['LopNr'].values)
    ys = list(dfy['LopNr'].values)

    sx  = set(xs)
    sy  = set(ys)
    cut = sx & sy
    ux  = sx - sy
    uy  = sy - sx

    print()
    # print('shape ' + namex + '\t\t', dfx.shape)
    # print('shape ' + namey + '\t\t', dfy.shape)
    # print('unique Lopnr ' + namex + '\t', len(xs))
    # print('unique Lopnr ' + namey + '\t', len(ys))

    print('common Lopnr\t\t\t', len(cut))
    print('Lopnr in ' + namex + ' only\t', len(ux))
    print('Lopnr in ' + namey + ' only\t', len(uy))
    print()

    ux = list(ux)
    uy = list(uy)
    ux.sort
    uy.sort
    return ux, uy


def readlopnr():
    dia = bio.read_generated_dia()
    per = bio.readperson()
    ctr = bio.readcontrol()
    inc = bio.readincare()
    nic = bio.readnicare()
    dru = bio.readdrug()
    tre = bio.readtreatment()
    sur = bio.readsurgery()
    cau = bio.readcause()

    data  = [dia, per, ctr, inc, nic, dru, tre, sur, cau]

    name = [
        'diagnos      ',
        'person       ',
        'kontrollgrupp',
        'sluten v_rd  ',
        '_ppen v_rd   ',
        'l_kemedel    ',
        'behandling   ',
        'kirurgi      ',
        'orsak        ',
    ]

    return data, name


def pairwise_lopnr_comparisions():

    data, name = readlopnr()

    for i in range(0, len(name)):
        for j in range(i+1, len(name)):
            print()
            print('--------------------------------------------------')
            print()
            print('\tComparing ' + name[i] + ' with ' + name[j])
            print()
            print('--------------------------------------------------')
            print()

            compare_lopnr(data[i], data[j], name[i], name[j])





""" -------------------------------
       
       4.  count amd list various types of diagnosis
           codes in care data
   
------------------------------------"""

"""
def is_icd10_class(x):
    if not common.isstr(x):
        return False
    if common.is_icd10(x):
        return False
    if len(x) < 3:
        return False
    if not x[0].isupper():
        return False
    return x[1].isdigit() and x[2].isdigit()
"""


def code_count(xs):
    if not isinstance(xs, str):
        return 0
    return len(xs.split())

def icd10_count(xs):
    if not isinstance(xs, str):
        return 0
    count = 0
    for x in xs.split():
        if common.is_icd10(x):
            # print(x)
            count += 1
    return count

def not_icd10_count(xs):
    if not isinstance(xs, str):
        return 0
    count = 0
    for x in xs.split():
        if not common.is_icd10(x):
            # print(x)
            count += 1
    return count

def icd10_class_count(xs):
    if not isinstance(xs, str):
        return 0
    count = 0
    for x in xs.split():
        if common.is_icd10_class(x):
            # print(x)
            count += 1
    return count

"""
def code_list(xs):
    if not isinstance(xs, str):
        return 0
    return len(xs.split())
"""

def count_and_print(df, table = False):
    dia  = 'DIAGNOS'
    dfc = copy.copy(df)
    dfc['code_count']        = df[dia].apply(code_count)
    dfc['icd10_count']       = df[dia].apply(icd10_count)
    dfc['not_icd10_count']   = df[dia].apply(not_icd10_count)
    dfc['icd10_class_count'] = df[dia].apply(icd10_class_count)
    nr_of_codes = dfc['code_count'].sum()
    nr_of_icd10 = dfc['icd10_count'].sum()
    nr_of_not_icd10 = dfc['not_icd10_count'].sum()
    nr_of_class_codes = dfc['icd10_class_count'].sum()

    if table:
        print('nr_of_lines\t', len(df))
        print('nr_of_codes\t', nr_of_codes)
        print('nr_of_icd10\t', nr_of_icd10)
        print('nr_of_not_icd10\t', nr_of_not_icd10)
        print('nr_of_icd10_class_codes\t', nr_of_class_codes)
    
    else:
            
    
        print('                nr_of_lines', len(df))
        print('                nr_of_codes', nr_of_codes)
        print('                nr_of_icd10', nr_of_icd10)
        print('            nr_of_not_icd10', nr_of_not_icd10)
        print('    nr_of_icd10_class_codes', nr_of_class_codes)


    """
    for c in df1[dia].values:
        print('\t', c)
    """


def print_dates(df, table = False):
    date = 'INDATUM'

    if table:

        print('first date\t', df[date].min())
        print('last date\t', df[date].max())

    else:

        print('           first date', df[date].min())
        print('            last date', df[date].max())
        

def icd10_class_list(xs):
    if not isinstance(xs, str):
        return []
    codes = []
    for x in xs.split():
        if common.is_icd10_class(x):
            codes += [x]
            #print(codes)
    return codes

def flat(xs):
    ys = []
    for x in xs:
        ys += x
    return ys

    

def print_class_codes(df):
    dia  = 'DIAGNOS'
    dfc = copy.copy(df)
    dfc['icd10_class'] = df[dia].apply(icd10_class_list)
    dfc['is_class'] = dfc['icd10_class'].apply(lambda x: x != [])
    dfc = dfc[dfc['is_class']]
    codes = np.unique(flat(list(dfc['icd10_class'].values)))
    for c in codes:
        print('\t', c)
 

def diagnosis_code_count(df, print_class = False, table = False):
    
    date = 'INDATUM'
    nr   = 'LopNr'
    icd10_start = np.datetime64('1998-01-01')

    """
    size0 = len(df)
    df = df.dropna().reset_index(drop=True)
    print('nr of empty lines:', size0- len(df))
    """
    
    df[date] = df[date].apply(bio.str2time)
    df = df.sort_values(date).dropna().reset_index(drop=True)

    df1 = df[df[date] < icd10_start]    
    df2 = df[df[date] >= icd10_start]

    print()    
    print('code counts before 1998_01_01:')
    print()
    
    print_dates(df1, table = table)
    count_and_print(df1, table = table)

    print()    
    print('code counts from 1998_01_01')
    print()
    
    print_dates(df2, table = table)
    count_and_print(df2, table = table)
    if print_class:
        print()
        print('    all icd10_class_codes:')
        print_class_codes(df2)

    print()
