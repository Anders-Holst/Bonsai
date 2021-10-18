#! /home/jan/anaconda3/bin/python


""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2013-03-26 version

    bonsai/src/lexicon.py

------------------------------------"""



import pandas as pd
import numpy as np
import bonsai_io as bio
import common
import copy


xcols    = ['ICD7', 'ICD9', 'text_C24_']
ycols    = ['ICD10', 'SNOMED']
lex_cols = xcols + ycols
tw_cols  = ['tw_icd10', 'tw_snomed']
dia_cols = ['LopNr'] + lex_cols




"""

   1   Common lexicon support for translation from
       ICD7, ICD9 to ICD10 and SNOMED

"""


def lookup(x, lex):
    for z in lex.values: 
        if (x[0] == z[0]) and (x[1] == z[1]) and common.isnull(x[2]) and common.isnull(z[2]):
            return z[3:5]
        if (x[0] == z[0]) and (x[1] == z[1]) and (x[2] == z[2]):
            return z[3:5]
    return [np.nan, np.nan]


def translate(df, lex):
    x = df[xcols]  # change dfrom df[1:4]
    y = lookup(x, lex)
    return y



"""

   2   Apply lexicon for translation from
       ICD7, ICD9 to ICD10 and SNOMED

"""


def apply_lex(df, lex, label):
    icd10  = 'icd10'+ label
    snomed = 'snomed'+ label
    df[icd10]  = df.apply(lambda d: translate(d, lex)[0], axis = 1)
    df[snomed] = df.apply(lambda d: translate(d, lex)[1], axis = 1)
    return df


def select_icd10(df, icd10):
    if common.isnull(df['ICD10']):
        return df[icd10]
    return df['ICD10']

def select_snomed(df, snomed):
    if common.isnull(df['SNOMED']):
        return df[snomed]
    return df['SNOMED']

"""

   3   Apply saved lexicon

"""



def fill_in_by_compiled_lex(df):
    lex = bio.read_code_lex()
    label = 'fill'    
    icd10  = 'icd10'+ label
    snomed = 'snomed'+ label
    df = apply_lex(df, lex, label)
    df['ICD10']  = df.apply(lambda d:select_icd10(d, icd10),  axis = 1)
    df['SNOMED'] = df.apply(lambda d:select_snomed(d, snomed), axis = 1)
    df = df.drop([icd10, snomed], axis=1)
    return df



"""

   4   Prune

"""

no_code = ['finns ej', 'finns ej som icd8kod', 'oklar kod']

def is_icd10(x, verb = False):
    if len(x) < 3:
        if verb: print(x, 'has less than 3 chars')
        return False
    if not x[0].isupper():
        if verb: print('\'' + x[0] + '\'' + ' is not an upper case letter')
        return False
    if not  x[1].isdigit():
        if verb: print('\'' + x[1] + '\'' + ' is not a digit')
        return False
    if not  x[2].isdigit():
        if verb: print('\'' + x[2] +'\'' + ' is not a digit')
        return False
    if len(x) > 6:
        if verb: print(x, 'has more than 6 chars')
        return False
    return True


def prune_icd10(c):
    if not isinstance(c, str) or c in no_code:
        return np.nan
    code = c
    if code == 'F':
        code = 'F00'
    if code == 'O':
        code = 'O50'
    if code == 'S':
        code = 'S83'
    code = code.replace('.', '')
    code = code.replace(' ', '')
    code = code[:3]
    # if not is_icd10(c, verb = True):
    if not is_icd10(c):
        print('\t' + c + '\t replaced by \t' +  code)
    return code

def prune_icd10_lex(df, col):
    df = df[df.columns[:2]]
    df.columns = [col, 'icd10']
    dfc = copy.copy(df)
    dfc['icd10'] = df['icd10'].apply(prune_icd10)
    dfc= dfc.sort_values(col).dropna().reset_index(drop=True)
    return dfc

def prune_icd10_lex_files(raw, lex, code):
    df = pd.read_csv(raw, sep=',', dtype = str)
    
    print('File:' + lex + '\n\tspecial icd10 replacements')
    print()
    df = prune_icd10_lex(df, code)
    print('saving file to:')
    print('\t', lex)
    df.to_csv(lex, index = False, sep = '\t')


def add_brackets(x):
    x = str(x)
    x = x.replace(' ', '')
    x = x.replace(',', '\',\'')
    s = '[\'' + x + '\']'
    return s


def add_brackets_to_file(f1, f2, col):
    df = pd.read_csv(f1, sep='\t', dtype = str)
    df[col] = df[col].apply(add_brackets)
    df.to_csv(f2, index = False, sep = '\t')

    
"""
def prune_icd10_lex_filesx(raw, lex, code):
    df = pd.read_csv(raw, sep=',', dtype = str)
    c1 = df.columns[1]
    df[c1] = df[c1].apply(str)
    print(df[0:3])
    df = df[df.columns[:2]]
    print(df[0:3])
"""


