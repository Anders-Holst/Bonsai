#! /usr/bin/env python3


""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE

    bonsai/src_v02/separate_tasks.py

    
------------------------------------"""

import numpy as np

from datetime import datetime as timeclass
from dateutil import parser
import pandas as pd

import bonsai_io as bio
import bonsai_time as btime
import common



""" -------------------------------

    finding zero dia group codes 
    
------------------------------------"""

def zero_grp():
    base = total.base_set()
    base = base[base[gs.names.dia_group] == '0']
    bio.save_tmp_file(base, 'base_grp_0')




""" -------------------------------

    mean diagnosis age
    
------------------------------------"""

def mean_diagnosis_age():
    da = 'Diagnos_lder'
    dia = bio.read_generated_dia()
    return dia[da].mean()


""" -------------------------------

    finding the data time spans 
    
------------------------------------"""


def spans():
    care_col  = 'INDATUM'
    dia_col   = 'DiagnosDat'
    cause_col = 'DODSDAT'
    drug_col  = 'EDATUM'

    inc = bio.readincare()
    t_inc = inc[care_col].dropna()
    t_inc = t_inc.sort_values()  # .apply(parser.parse)

    nic = bio.readnicare()
    t_nic = nic[care_col].dropna()
    t_nic = t_nic.sort_values() # .apply(parser.parse)

    dia = bio.read_generated_dia()
    t_dia = dia[dia_col].dropna() # .apply(parser.parse)

    oc  = bio.readoutcome()
    t_oc = oc[dia_col].dropna() # .apply(parser.parse)

    oc = oc[oc[cause_col] != 0]
    t_life = oc[cause_col].dropna().apply(str).apply(parser.parse).apply(common.time2str)
    # t_life = t_life[t_life[cause_col] != '0']
    
    drug = bio.readdrug()
    t_drug = drug[drug_col].dropna()# .apply(parser.parse).apply(common.time2str)


    print()
    print('incare span: ', t_inc.min(), '\t', t_inc.max())
    print('nicare span: ', t_nic.min(), '\t', t_nic.max())
    print('dia span:    ', t_dia.min(), '\t', t_dia.max())
    print('oc span:     ', t_oc.min(), '\t', t_oc.max())
    print('life span:   ', t_life.min(), '\t', t_life.max())
    print('drug span:   ', t_drug.min(), '\t', t_drug.max())
    print()

    #t = list(inc[care_col].values)


""" -------------------------------

    generation of diagnosis codes lexicon for incare  
    
------------------------------------"""

xdir = '/home/jan/bonsai/data/extracted/'
ks68_file  = xdir + 'diagnoskoder_20190425_sv_70_86.csv'
icd9_file  = xdir + 'diagnoskoder_20190425_sv_86_97.csv'


def heads_list(c):
    nc = []
    for ci in c:
        if common.isstr(ci):
            nc += [ci[0]]
    return nc
    
def icd10_heads_list(df):
    c1 = df['ICD10_1']
    c2 = df['ICD10_2']
    c  = [c1] + [c2]
    return heads_list(c)

def gen_incare_lex():
    c1 = 'ICD10_1'
    c2 = 'ICD10_2'
    cols = ['code', c1, c2]
    f1  = ks68_file
    f2  = icd9_file
    df1 = pd.read_csv(f1, dtype = str, sep='\t')
    df2 = pd.read_csv(f2, dtype = str, sep='\t')
    df2 = df2[ ~df2[c1].isnull() | ~df2[c2].isnull() ]
    df1.columns = cols
    df2.columns = cols

    df1['ICD10'] = df1.apply(icd10_heads_list, axis = 1)
    df2['ICD10'] = df2.apply(icd10_heads_list, axis = 1)
    df = pd.concat([df1, df2])
    df = common.rmcols(df, [c1, c2])

    bio.save_local_input(df, 'incare_lex')



""" -------------------------------

    generation of diagnosis codes lexicon for causes   
    
------------------------------------"""


xdir = '/home/jan/bonsai/data/extracted/'
cause_file = xdir + 'diagnoskoder_20190425_orsak.csv'


def gen_cause_lex():

    c1 = 'ICD10_1'
    c2 = 'ICD10_2'
    cols = ['code', c1, c2]

    f = cause_file
    df = pd.read_csv(f, dtype = str, sep=',')
    df = df[ ~df[c1].isnull() | ~df[c2].isnull() ]
    df.columns = cols
    df['ICD10'] = df.apply(icd10_heads_list, axis = 1)
    df = common.rmcols(df, [c1, c2])

    df8 = df[ df['code'].str.contains('#')]
    df9 = df[ ~df['code'].str.contains('#')]

    df8.is_copy = False
    df8['code'] = df8['code'].str[1:] 
    
    bio.save_local_input(df8, 'cause_icd8_lex')
    bio.save_local_input(df9, 'cause_icd9_lex')
