#! /home/jan/anaconda3/bin/python


""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2013-03-26 version

    bonsai/src_v02/bonsai_time.py
    
    support for managing of time data 
    and time analysis


------------------------------------"""



from datetime import datetime as timeclass
from datetime import timedelta as delta
from dateutil import parser
import pandas as pd
import numpy as np

import common
import global_settings as gs
import bonsai_io as bio


def head(xs, lex8, lex9, n = 1):
    ys = []
    for x in xs:
        if isinstance(x, str):
            cx = common.lookup_icd10(x, lex8, lex9)
            for c in cx:
                ys +=  [c[:n]]
    return ys


def split_and_head(ys, recidiv_icd10_class, lex8, lex9, n = 1):
    
    cs = []
    for y in ys.values:
        if not common.isnull(y):
            cs = np.append(cs, y.split())
    cs = list(cs)
    ys = np.unique(head(cs, lex8, lex9, n = n))
    ys = list(set(ys) - set([recidiv_icd10_class]))
    return ys

# here we can instead return the difference

def after(x, y, d = 0):
    if common.notnull(x) and common.notnull(y):
        return parser.parse(x) > parser.parse(y) +  delta(days = d)
    return False

# here we can instead return the time and difference

def times_after(sz, y, d = 0):
    L = []
    for z in sz:
        if after(z, y, d):
            L += [z]
    return L

    
def first_time(zs):
    L = list(map(common.str2time, zs))
    if L == []:
        return '0'
    # t = min(list(map(common.str2time, zs)))
    return common.time2str(min(L))


def first_after(zs, y, d = 0):
    return first_time(times_after(zs, y, d))


def first_time_compression(df, xcol, ycol, zcol, fcol, lex8, lex9, gap = 1826, split = True, n = 1):

    # xcol = 'LopNr'
    # ycol = 'DIAGNOS'
    # fcol = 'first_incare'

    data = []
    df = df.dropna()
    if df.empty:
        return df

    xs =  df[xcol].drop_duplicates()

    print()
    print('first_time_compression')
    print()
    i = 0

    print('nr of LopNr:s',  len(xs))
    for x in xs:
        i += 1
        if (i % 100) == 0:
            print(i)
        
        dx = df[ df[xcol] == x ]  
        diagnos_dat   = dx['DiagnosDat'].values[0]
        recidiv_icd10 = dx['ICD10'].values[0]    
        recidiv_icd10_class = recidiv_icd10[:n]
        sz = dx[zcol].drop_duplicates()
        yz_list = []

        for z in sz:
            dz = dx[ dx[zcol] == z ]  
            ys = dz[ycol].drop_duplicates()


            if split:
                ys = split_and_head(ys, recidiv_icd10_class, lex8, lex9, n = n)
                
            for y in ys:
                yz_list += [[y, z]]

        if not (yz_list == []):
            dyz = pd.DataFrame(yz_list)
            dyz.columns = [ycol, zcol]      
            sy = dyz[ycol].drop_duplicates()
            yminz_list = []
            for y in sy:
                dy = dyz[ dyz[ycol] == y ]  
                times = dy[zcol].values
                z = first_after(times, diagnos_dat, gap)
                if z != '0':
                    yminz_list += [[y, z]]
            data = data + [ [x] + [yminz_list] ]

    dz = pd.DataFrame(data)
    if not dz.empty:
        dz.columns = [xcol, fcol]
    return dz

def first_time_aho(df, base, xcol, ycol, zcol, fcol, lex8, lex9, gap = 1826, split = True, n = 1):

    # xcol = 'LopNr'
    # ycol = 'DIAGNOS'
    # fcol = 'first_incare'

    data = []

    xs =  base[xcol]

    print()
    print('first_time_compression aho')
    print()
    i = 0

    df.set_index(xcol,drop=False,inplace=True)

    print('nr of LopNr:s',  len(xs))
    for x in xs:
        i += 1
        if (i % 100) == 0:
            print(i)
        
        if x not in df.index:
            continue

        dx = df.loc[x]  
        if type(dx) != pd.DataFrame:
            dx = pd.DataFrame([dx])
        diagnos_dat   = base['DiagnosDat'][x]
        recidiv_icd10 = base['ICD10'][x]
        recidiv_icd10_class = recidiv_icd10[:n] if isinstance(recidiv_icd10, str) else ""

        yz_list = []
        for j in range(len(dx)):
            yz_list += split_and_head_1(dx.iloc[j], ycol, zcol, recidiv_icd10_class, lex8, lex9, n = n)
                
        if not (yz_list == []):
            dyz = pd.DataFrame(yz_list)
            dyz.columns = [ycol, zcol]      
            sy = dyz[ycol].drop_duplicates()
            yminz_list = []
            for y in sy:
                dy = dyz[ dyz[ycol] == y ]  
                times = dy[zcol].values
                z = first_after(times, diagnos_dat, gap)
                if z != '0':
                    yminz_list += [[y, z]]
            data = data + [ [x] + [yminz_list] ]

    dz = pd.DataFrame(data)
    if not dz.empty:
        dz.columns = [xcol, fcol]
    return dz



def all_times(df, xcol, ycol, zcol, acol, lex8, lex9, split = True, n = 1):
    
    # xcol = 'LopNr'
    # ycol = 'DIAGNOS'
    # acol = 'all_incare'
    
    data = []
    df = df.dropna()
    if df.empty:
        return df
    
    xs =  df[xcol].drop_duplicates()

    print()
    print('all_time_compression')
    print()
    i = 0
    
    for x in xs:

        i += 1
        if (i % 100) == 0:
            print(i)

        dx = df[ df[xcol] == x ]  
        diagnos_dat   = dx['DiagnosDat'].values[0]
        recidiv_icd10 = dx['ICD10'].values[0]
        recidiv_icd10_class = recidiv_icd10[:n]
        sz = dx[zcol].drop_duplicates()
        yz_list = []

        for z in sz:
            dz = dx[ dx[zcol] == z ]  
            ys = dz[ycol].drop_duplicates()
            
            if split:
                ys = split_and_head(ys, recidiv_icd10_class, lex8, lex9, n = n)
            
            for y in ys:
                yz_list += [[y, z]]
                
        if not (yz_list == []):
            dyz = pd.DataFrame(yz_list)
            dyz.columns = [ycol, zcol]      
            sy = dyz[ycol].drop_duplicates()

            yallz_list = []
            for y in sy:
                dy = dyz[ dyz[ycol] == y ]  
                times = list(dy[zcol].values)
                yallz_list += [[y, times]]
            data = data + [ [x] + [yallz_list] ]

    dz = pd.DataFrame(data)
    if not dz.empty:
        dz.columns = [xcol, acol]    
    return dz

def split_and_head_1(row, ycol, zcol, recidiv_icd10_class, lex8, lex9, n = 1):
    if common.isnull(row[ycol]) or common.isnull(row[zcol]):
        return []
    cs = row[ycol].split()
    ys = np.unique(head(cs, lex8, lex9, n = n))
    ys = list(set(ys) - set([recidiv_icd10_class]))
    return [[y, row[zcol]] for y in ys]

def all_times_aho(df, base, xcol, ycol, zcol, acol, lex8, lex9, split = True, n = 1):
    
    # xcol = 'LopNr'
    # ycol = 'DIAGNOS'
    # acol = 'all_incare'
    
    data = []
    
    xs =  base[xcol]

    print()
    print('all_time_compression aho')
    print()
    i = 0

    df.set_index(xcol,drop=False,inplace=True)

    for x in xs:

        i += 1
        if (i % 100) == 0:
            print(i)

        if x not in df.index:
            continue

        dx = df.loc[x]  
        if type(dx) != pd.DataFrame:
            dx = pd.DataFrame([dx])
        diagnos_dat   = base['DiagnosDat'][x]
        recidiv_icd10 = base['ICD10'][x]
        recidiv_icd10_class = recidiv_icd10[:n] if isinstance(recidiv_icd10, str) else ""

        yz_list = []
        for j in range(len(dx)):
            yz_list += split_and_head_1(dx.iloc[j], ycol, zcol, recidiv_icd10_class, lex8, lex9, n = n)
                
        if not (yz_list == []):
            dyz = pd.DataFrame(yz_list)
            dyz.columns = [ycol, zcol]      
            sy = dyz[ycol].drop_duplicates()

            yallz_list = []
            for y in sy:
                dy = dyz[ dyz[ycol] == y ]  
                times = list(dy[zcol].values)
                yallz_list += [[y, times]]
            data = data + [ [x] + [yallz_list] ]

    dz = pd.DataFrame(data)
    if not dz.empty:
        dz.columns = [xcol, acol]    
    return dz
