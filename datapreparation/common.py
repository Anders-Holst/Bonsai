#! /usr/bin/env python3


""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2019-08-29 version

    bonsai/src/common.py
    common support 

------------------------------------"""

import pandas as pd
import numpy as np
from datetime import datetime as timeclass
from datetime import timedelta as delta
from dateutil import parser
import ast
import copy

from scipy.stats import beta
# from scipy.stats import expon
# from scipy.stats import gamma


""" ---------------------------

   boolean functions

------------------------- """


def isstr(x):
    return isinstance(x, str)


def isnull(a):
    return not(isinstance(a, str))


def notnull(x):
    if isstr(x):
        return (x != '0')
    return False
 

def isarray(x):
    return (isinstance(x, list) or isinstance(x, np.ndarray))

# what is an icd 10 code, we require that it starts
# with an uppercase letter followed by three numbers


def is_icd10(x):
    if not isstr(x):
        return False
    if len(x) < 4:
        return False
    if not x[0].isupper():
        return False
    return x[1].isdigit() and x[2].isdigit() and x[3].isdigit()




""" ---------------------------

   conversions

------------------------- """



def str2number(x, dec = 5):
    s = str(x)
    if  (',' in s):
        x = float(s.replace(',','.'))
    return round(float(x), dec)


def str2list(s):
    return ast.literal_eval(s)

def bool2str(x):
    if isinstance(x, bool) or isinstance(x, int):
        return str(int(x))
    return x

def df_bool2str(df, cols):
    for c in cols:
        df[c] = df[c].apply(bool2str)
    return df




""" ---------------------------

   lists

   [[['JF', 'PJ'], 41], [['JE'], 0]] to
   [['JF', 41], ['PJ', 41], ['JE', 0]]


------------------------- """

def code_compress(u):
    ru = []
    for x in u:
        rx = []
        [cs, d] = x
        for c in cs:
            rx += [[c, d]]
        ru += rx
    return ru


""" ---------------------------

   pandas

------------------------- """

def rmcol(df, c):
    return df.drop([c], axis=1)

def rmcols(df, cs):
    for c in cs:
        df = df.drop([c], axis=1)
    return df


def change_col_names(df, old_name, new_name):
    df.columns = [new_name if x==old_name else x for x in df.columns]
    return df


def to_each_col_apply(df, cols, f):
    for c in cols:
        df[c] = df[c].apply(f)



def compression(df, xcol, vcol):
    data = []
    xs =  df[xcol].drop_duplicates()
    for x in xs:
        dx = df[ df[xcol] == x ]  
        codes = []
        for line in dx[vcol].values:
            if isinstance(line, list):
                codes += line
        codes = np.unique(codes)
        data += [[x, codes]]
        
    d = pd.DataFrame(data)
    d.columns = [xcol, vcol]
    d[vcol] = d[vcol].apply(list)
    return d



def compression_with_dates(df, xcol, vcol, dcol):
    data = []
    xs =  df[xcol].drop_duplicates()
    for x in xs:

        dx = df[ df[xcol] == x ]
        L0 = len(dx)
        dx = dx[dx[dcol] == dx[dcol]] # remove NaT values
        L1 = len(dx)
        if L0 != L1:
            print(L0-L1, 'non-valid surgery dates for this LopNr')
            print(dx)
            print()
 
        codes = []
        dxcopy = copy.copy(dx) 
        dxcopy['Datum'] = dx['Datum'].apply(lambda x: x.strftime('%Y-%m-%d'))
        dx = dxcopy
        for values_and_date in dx[[vcol, dcol]].values:
            values_and_date = [values_and_date[0]] + [values_and_date[1]]

            if isinstance(values_and_date, list):
                codes += [values_and_date]
        data += [[x, codes]]
        
    d = pd.DataFrame(data)
    d.columns = [xcol, vcol]
    d[vcol] = d[vcol].apply(list)
    return d




def compression_date_limit(df, xcol, vcol, dcol, date_dict):
    data = []
    xs =  df[xcol].drop_duplicates()
    for x in xs:
        dx = df[ df[xcol] == x and df[dcol] <= date_dict[x]]  
        codes = []
        for line in dx[vcol].values:
            if isinstance(line, list):
                codes += line
        codes = np.unique(codes)
        data += [[x, codes]]
        
    d = pd.DataFrame(data)
    d.columns = [xcol, vcol]
    d[vcol] = d[vcol].apply(list)
    return d





""" ---------------------------

   lexicon (same code as in merge.py)

------------------------- """



def look_up_entry(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        return dfe[value_col].values[0] 
    return str(0)


def look_up_value(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        return dfe[value_col].values[0] 
    return np.nan


def look_up_list_value(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        return str2list( dfe[value_col].values[0] )  
    return str(0)


def look_up_codes_1(x, lex, col):
    keys = set(lex[col].values)
    if x in keys:
        return look_up_list_value(x, lex, col, 'icd10')
    elif x[0].isalpha() and x[1:] in keys:
        return look_up_list_value(x[1:], lex, col, 'icd10')
    else:
        return []

def look_up_codes(x, lex, col):
    keys = set(lex[col].values)
    if x in keys:
        return look_up_list_value(x, lex, col, 'icd10')
    return []


def look_up_list(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if not dfe.empty:
        return dfe[value_col].values[0]
    return []


#def lookup_icd10(x, lex8, lex9):
#    lex8.columns = ['code', 'ICD10']
#    lex9.columns = ['code', 'ICD10']
#    lex = pd.concat( [lex8, lex9])
#    return look_up_codes(x, lex)

def lookup_icd10(x, lex8, lex9):
    res = look_up_codes(x, lex8, 'icd8') if lex8 is not False else []
    if res:
        return res
    res = look_up_codes(x, lex9, 'icd9') if lex9 is not False else []
    if res:
        return res
    if is_icd10(x) or is_icd10_class(x):
        return [x]
    else:
        return []

def add_list_cols(df1, df2, xc, ycs):

    # print(df1[:5]) # zzz
    # print(df2[:5]) # zzz
    
    for yc in ycs:
        df1[yc] = df1[xc].apply(lambda x:look_up_list(x, df2, xc, yc))
    return df1


def look_up_all_values(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if dfe.empty:
        return []
    return list(dfe[value_col].values)

"""
def look_up_all_values_and_dates(entry, df, entry_col, value_col, date_col):
    dfe = df[df[entry_col] == entry]
    return dfe[[value_col, date_col]]
"""


def look_up_all_values_date_limit(entry, df, entry_col, value_col, date_col, date_dict):
    dfe = df[df[entry_col] == entry and df[date_col] <= date_dict[entry]]
    if dfe.empty:
        return []
    return list(dfe[value_col].values)


def look_up_nonnan(entry, df, entry_col, value_col):
    dfe = df[df[entry_col] == entry]
    if dfe.empty:
        return []
    v = list( set(dfe[value_col].values) - set([np.nan]) )
    return v


""" ---------------------------

   Time

   t: timeclass
   x: time as a string (str)

   Here we use the string '0' as an undefined value

------------------------- """

def is_time(t):
    # return isinstance(t, datetime.datetime)
    return isinstance(t, timeclass) 

def str2time(x):
    if notnull(x):
        return parser.parse(x)
    return x

"""
def full_year2time(x):
    return str2time(x + '-06-15')
"""

def time2str(t):
    if is_time(t):
        return t.strftime("%Y-%m-%d")
    return t


def time_diff(t1, t0):
    if is_time(t1) and is_time(t0):
        return (t1-t0).total_seconds()
    return np.nan  


def time_diff_in_seconds(x, y):
    if notnull(x) and notnull(y):
        if x[-2:]=='00':
            x = x[:-2] + "15"
        if y[-2:]=='00':
            y = y[:-2] + "15"
        tx = parser.parse(x)
        ty = parser.parse(y)
        return (tx-ty).total_seconds()
    return np.nan   

def time_diff_in_days(x, y):
    s = time_diff_in_seconds(x, y)
    if np.isnan(s) or (s == '0'):
        return s
    return round( s/ (60*60*24) )


def time_diff_in_years(x, y):
    s = time_diff_in_days(x, y)
    if np.isnan(s) or (s == '0'):
        return s
    return round( s/ 365.24, 2 )

def time_min(x, y):
    return time2str(min(parser.parse(x), parser.parse(y)))

def time_max(x, y):
    return time2str(max(parser.parse(x), parser.parse(y)))

def add_days(x, nr_of_days):
    # print(time2str(parser.parse(x)))
    return time2str(parser.parse(x) +  delta(days = nr_of_days))

def add_years(x, nr_of_years):
    nr_of_days = 365.24 * nr_of_years
    return add_days(x, nr_of_days)

def add_years_to_full_year(x, nr_of_years):
    x = x + '-07-02'
    return add_years(x, nr_of_years)


""" ---------------------------

   two sample test

------------------------- """


def exp2sample_prob(nx, ny, sx, sy):
    if sx == 0 or sy == 0 or nx == 0 or ny == 0:
        return 0.5
    h = sy / (sx + sy)
    return beta.cdf(h, ny, nx)


def exp2sample_logprob(nx, ny, sx, sy):
    if sx == 0 or sy == 0 or nx == 0 or ny == 0:
        return 0.5
    h = sy / (sx + sy)
    return beta.logcdf(h, ny, nx)


def exp2sample_logsf(nx, ny, sx, sy):
    if sx == 0 or sy == 0 or nx == 0 or ny == 0:
        return 0.5
    h = sy / (sx + sy)
    return beta.logsf(h, ny, nx)


""" ---------------------------

   fill in icd10 class codes

------------------------- """




def is_icd10_class(x):
    if not isstr(x):
        return False
    if is_icd10(x):
        return False
    if len(x) < 3:
        return False
    if not x[0].isupper():
        return False
    return x[1].isdigit() and x[2].isdigit()

def icd10_fill(xs, s = '0*'):
    if not isinstance(xs, str):
        return xs
    ys = ''
    for x in xs.split():
        if is_icd10_class(x):
            ys += ' ' + x[:3] + s + x[3:]
            # print('\t', x, '\t', x[:3] + s + x[3:])
        else:
            ys += ' ' + x 
    return ys[1:]
