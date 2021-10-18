#! /usr/bin/env python3


""" -------------------------------

    bonsai_io.py

    Copyright (C) 2019 RISE
    This code was produced by RISE
    The 2019-08-29 version

    bonsai/src/bonsai_io.py

------------------------------------"""


import copy
import pandas as pd
import numpy as np
from datetime import datetime as time
from datetime import timedelta as delta
from dateutil import parser


import global_settings as gs
import common
import diagnose


""" ----------------------------------------

   common support

---------------------------------------"""



def isstr(x):
    return isinstance(x, str)

def str2number(x):
    s = str(x)
    if  (',' in s):
        x = float(s.replace(',','.'))
    return round(float(x), 5)

def str2time(x):
    # print(x)
    if isstr(x):
        return parser.parse(x)
    return x

def save_generated_file(df, name):
    f = gs.places.generated_out(name)
    print('saving ' + name + ' file to:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')

def read_generated_file(name):
    f = gs.places.generated_out(name)
    print('reading file:', f)
    df = pd.read_csv(f, sep='\t')
    return df
    

def save_tmp_file(df, name):
    f = gs.places.tmp_out(name)
    print('saving ' + name + ' file to:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')

def save_local_tmp_file(df, name):
    f = gs.places.local_tmp_out(name)
    print('saving ' + name + ' file to:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')


def save_local_result(df, name):
    f = gs.places.local_result(name)
    print('saving ' + name + ' file to:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')


def save_local_input(df, name):
    f = gs.places.local_input + name + '.csv'
    print('saving ' + name + ' file to:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')



""" ----------------------------------------

    lexicon 

---------------------------------------"""

def read_code_lex():
    lex = pd.read_csv(gs.places.code_lex, dtype = str, sep='\t')
    return lex


def read_incare_lex():
    name = 'incare_lex'
    f = gs.places.local_input + name + '.csv'
    print('reading file:', f)
    lex = pd.read_csv(f, dtype = str, sep='\t')
    return lex


def read_care_icd8_lex():
    f = gs.places.care_icd8_lex
    print('reading file:', f)
    lex = pd.read_csv(f, dtype = str, sep='\t')
    return lex


def read_care_icd9_lex():
    f = gs.places.care_icd9_lex
    print('reading file:', f)
    lex = pd.read_csv(f, dtype = str, sep='\t')
    return lex

def read_cause_icd8_lex():
    f = gs.places.cause_icd8_lex
    print('reading file:', f)
    lex = pd.read_csv(f, dtype = str, sep='\t')
    return lex


def read_cause_icd9_lex():
    f = gs.places.cause_icd9_lex
    print('reading file:', f)
    lex = pd.read_csv(f, dtype = str, sep='\t')
    return lex

""" ----------------------------------------

    dia

    reading the diagnosis file, the original 
    and saved filled in versions

---------------------------------------"""


def read_all_original_dia():
    """
    reading the the non-null LopNr rows from the
    original file helena_lev_Diagnos_v2.txt
    modified only by removing a ',' sign
    """

    f  = gs.places.diagnose_m
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['LopNr'].isnull() ]
    # df = df[gs.places.diagnose_selection]
    df = df.sort_values('LopNr').reset_index(drop=True)
    return df

def read_original_dia():
    """
    reading the non-null LopNr rows and the selected 
    columns from original file helena_lev_Diagnos_v2.txt 
    modified only by removing a ',' sign
    """

    f  = gs.places.diagnose_m
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.diagnose_selection]
    df = df.sort_values('LopNr').reset_index(drop=True)
    return df


def save_generated_dia(dia):
    print('saving generated diagnostics file to:')
    print('\t', gs.places.generated_dia)
    dia.to_csv(gs.places.generated_dia, index = False, sep = '\t')

"""
def parse(s):
    print(s, type(s))
    print(str(s))
    return parser.parse(str(s))
"""

def read_generated_dia():

    f = gs.places.generated_dia
    print('reading file:', f)
    dia = pd.read_csv(f, dtype = str, sep='\t')
    da = 'Diagnos_lder'
    dia[da] = dia[da].apply(common.str2number)
    """
    for col in gs.places.diagnose_time_cols:
        dia[col] = dia[col].apply(str2time)
    """
    return dia

def save_tmp_dia(dia):
    print('saving generated diagnostics groups file to:')
    print('\t', gs.places.tmp_dia)
    dia.to_csv(gs.places.tmp_dia, index = False, sep = '\t')



""" ----------------------------------------

    dia groups 

    reading and saving versions of the 
    diagnosis groups file

---------------------------------------"""

def group_nr(s):
    return int(s.split()[1][0:-1])

def rm_dot(s):
    return s.replace('.','')
 
def read_orig_diagroups():
    grp = gs.names.input_data_group
    df = pd.read_csv(gs.places.orig_dia_groups, dtype = str, sep='\t')
    cols = list(df.columns[4:6]) + list([df.columns[1]]) 
    df = df[cols]
    df.columns = ['ICD10',  'SNOMED', grp]
    df[grp] = df[grp].apply(group_nr)
    df['ICD10'] = df['ICD10'].apply(rm_dot)
    return df
   
def save_diagroups():
    df = read_orig_diagroups()
    df = df.sort_values(['ICD10', 'SNOMED']).reset_index(drop=True)
    print('saving generated diagnostics groups file to:')
    print('\t', gs.places.dia_groups)
    df.to_csv(gs.places.dia_groups, index = False, sep = '\t')


def read_diagroups():
    f = gs.places.dia_groups
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    return df




""" ----------------------------------------

    person 
    reading the person data file

---------------------------------------"""


def k_n(x):
    if not(x == 'flicka') and not(x == 'pojke'):
        return 'sex unknown'
    return x

def vital(x):
    return {
        'levande':'alive',
        'Död':'dead',
        'utvandrade':'emigrated',
        'okänt':'emigrated'
    }[x]



def readperson():
    time_cols = gs.places.person_time_cols
    cols = gs.places.person_cols
    
    f  = gs.places.person
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.person_selection]    
    for c in time_cols:
        df[c] = df[c].apply(str2time)  # may destroy the values of other cols
    df['K_n'] = df['K_n'].apply(k_n)
    
    # df['last_date'] = df['Registerdata_inh_mtade__CR_Bef__']
    df['VitalStatus'] = df['VitalStatus'].apply(vital)
    
    df.columns = cols
    return df



""" ----------------------------------------

    drug
    this concerns the outcome, drugs for everyone

---------------------------------------"""


def readdrug(n = 1):
    f  = gs.places.drug
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    # df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[ df['LopNr'].str.contains('\d', regex=True) ]
    df = df[gs.places.drug_selection]  
    df['ATC'] = df['ATC'].str[:n].apply(str) # zzz
    df = df.drop_duplicates()
    return df


""" ----------------------------------------

    cytostatica

---------------------------------------"""


def readcytocodes():
    f  = gs.places.cytocodes
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['BehNr'].isnull() ]
    df = df[gs.places.cyto_selection]
    df = df.drop_duplicates()
    return df
  

def readdrugclasses():
    f  = gs.places.drug_classes
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[gs.places.drug_classes_selection]  
    df = df.drop_duplicates()
    return df



""" ----------------------------------------

    treatement
    reading the tratement file

    not complete yet 2019-04-02

---------------------------------------"""


def readtreatment():
    f  = gs.places.treatment
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t')

  
    
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.treatment_selection]

    
    # df = df.dropna() # will drop all but 18 LopNr:s


    df = common.change_col_names(df, 'StartDatum_Cyt', 'cyto_date')
    df = common.change_col_names(df, 'StartDatum_Stamcell', 'stem_date')
    df = common.change_col_names(df, 'KirurDat', 'surg_date')
    df = common.change_col_names(df, 'StartDatum_Radioterapi', 'radio_date')
    df = df.drop_duplicates().reset_index(drop=True)

    # xs = np.unique(df['LopNr'].values)
    # print('nr of LopNr', len(xs))
  
    
    return df





""" ----------------------------------------

    incare and nicare
    institutional and non-institutional care
    reading the files

---------------------------------------"""

def readorigincare():
    f  = gs.places.orig_incare
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.incare_selection]  
    return df


def readorignicare():
    f  = gs.places.orig_nicare
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.nicare_selection]  
    df = df.drop_duplicates() 
    return df



def readincare():
    f  = gs.places.incare
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    return df


def readnicare():
    f  = gs.places.nicare
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    return df



def modify_care(df):
    
    date = 'INDATUM'
    dia = 'DIAGNOS'
    icd10_start = np.datetime64('1998-01-01')
    df[date] = df[date].apply(str2time)
    df = df.sort_values(date).dropna().reset_index(drop=True)    
    df1 = df[df[date] < icd10_start]    
    df2 = df[df[date] >= icd10_start]
    df2c = copy.copy(df2)
    df2c[dia] = df2[dia].apply(common.icd10_fill)
    df = pd.concat( [df1, df2c])
    return df

def modify_incare():
    df = readorigincare()
    df = modify_care(df)
    f  = gs.places.incare
    print('saving to file:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')
    

def modify_nicare():
    df = readorignicare()
    df = modify_care(df)
    f  = gs.places.nicare
    print('saving to file:')
    print('\t', f)
    df.to_csv(f, index = False, sep = '\t')

""" ----------------------------------------

    control group

---------------------------------------"""

def readcontrol():
    f  = gs.places.control
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    df = df[gs.places.control_selection]  
    df.columns = gs.places.control_column_names

    # fc = 'FoddAr'
    # df[fc] = df[fc].apply(parser.parse)
    # df = df[ ~df['LopNr'].isnull() ]
    # df = df[gs.places.control_selection]  
    return df


""" ----------------------------------------

    cause

    because of sever problems with pandas: a datetime object 
    in one column interfers with operations on other columns,
    the readcause function manages the cause list one letter 
    construction

---------------------------------------"""



def old_add_to_cause_list(d, c, lex8, lex9):
    recidiv_icd10 = d['ICD10'] # do no add first letter of recidiv

    icd = d['ICD']
    L = d['life_events']
    x = d[c]
    if L == 'null': #  L is null from start 
        L = []

    if not isinstance(x, str):
        return L

    if x == recidiv_icd10:
        return L

    if icd == '10':    
        L = list(np.unique( L + [x[0]] ))
        return L

    if icd == '9':
        cx = common.look_up_codes(x, lex9)
    if icd == '8':
        cx = common.look_up_codes(x, lex8)

    for c in cx:
        L +=  [c[0]]
    L = list(np.unique( L ))
    return L


def add_class(L, x, rc, nr = 3):
    xc = x[:nr]
    if not xc == rc:
        L = list(np.unique( L + [xc] ))
    return L


def add_to_cause_list(d, c, lex8, lex9, nr = 3):
    # nr is the nr of icd10 chars of concern
    
    recidiv_icd10 = d['ICD10'] # do no add recidiv class diagnosis
    recidive_class = recidiv_icd10[:nr]

    icd = d['ICD']
    L = d['life_events']
    x = d[c]
    
    if L == 'null': #  L is null from start 
        L = []

    if not isinstance(x, str):
        return L
    
    if icd == '10':
        if common.is_icd10(x):
            L = add_class(L, x, recidive_class, nr = nr)
        elif not x[0].isalpha() and not x[0].isdigit() and common.is_icd10(x[1:]):
            L = add_class(L, x[1:], recidive_class, nr = nr)
        return L

    if icd == '9':
        cx = common.look_up_codes_1(x, lex9, 'icd9')
        
    if icd == '8':
        cx = common.look_up_codes_1(x, lex8, 'icd8')


    if common.isarray(cx):
        for c in cx:
            L = add_class(L, c, recidive_class, nr = nr)
        
    elif common.isstr(cx):
        L = add_class(L, cx, recidive_class, nr = nr)
    
    return L
  

def readcause(df_icd10):
    
    # read cause and add a column 'life_events' 
    # with a list of first letter of non-recidiv causes
    
    cfile  = gs.places.cause
    print('reading file:', cfile)
    df = pd.read_csv(cfile, dtype = str, sep='\t', encoding='latin1')
    df = df[gs.places.cause_selection]  
    df = diagnose.add_cols_or_zero(df, df_icd10, 'LopNr', ['ICD10'])
    df['life_events'] = 'null'

    lex8 = read_cause_icd8_lex()
    lex9 = read_cause_icd9_lex()

    for c in gs.places.morsak:
        f = lambda d:add_to_cause_list(d, c, lex8, lex9)  # use cause lex here
        df['life_events'] = df.apply(f, axis = 1)
    cols =  ['LopNr'] + gs.places.cause_added_cols
    df = df[cols] 
    """
    this destroys it all 


    !!
    for col in gs.places.cause_time_cols:
        df[col] = df[col].apply(str2time)
    """

    return df
 

""" ----------------------------------------

    surgery

---------------------------------------"""

def readsurgery():
    f  = gs.places.surgery
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', 
                     keep_default_na = False,
                     na_values = [''])


    df = df[ ~df['LopNr'].isnull() ]
    # remove the non-number LopNr in the surgery file
    df = df[ df['LopNr'].str.contains('\d', regex=True) ] 
    df = df[gs.places.surgery_selection] 
    df = df.drop_duplicates()
    df.columns = gs.places.surgery_column_names
    df['Datum'] = df['Datum'].apply(str2time)
    return df
 
def read_generated_surgery():
    df = read_generated_file('surgery')
    df['LopNr'] = df['LopNr'].apply(str)
    df['surgery'] = df['surgery'].apply(common.str2list)
    return df


""" ----------------------------------------

    stralung

---------------------------------------"""


def readradiocodes():      # BehNr > Target_Kod
    f = gs.places.radio
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    df = df[['LopNr', 'BehNr', 'Target_Kod']]
    return df

def readradioclasses():    # Target_kod > Kategorier
    f = gs.places.radio_classes
    print('reading file:', f)
    df = pd.read_csv(f, dtype = str, sep='\t', encoding='latin1')
    df = common.rmcol(df, 'Target')

    # df = df[['Target_Kod', 'k1', 'k2', 'k3', 'k4']]
    return df




""" -------------------------------

    special missions

    reading the  columns Target_I and Target_II
    otherwise not used, from the treatement file 

------------------------------------"""

def readtarget():
    f  = gs.places.treatment
    df = pd.read_csv(f, dtype = str, sep='\t')
    df = df[ ~df['LopNr'].isnull() ]
    df = df[gs.places.treatment_tagets]  
    return df




""" -------------------------------

    evolving reads

------------------------------------"""


# use bio.read_generated_file(gs.names.outcome_file) instead

def readoutcome():
    col = 'first_incare'
    name = gs.names.outcome_file
    df   = read_generated_file(name)
    df[col] = df[col].apply(common.str2list)
    return df


def read_final(name):
    ei_col = 'event_time_inc'
    ci_col = 'event_time_death'
    df = read_generated_file(name)
    df['LopNr'] = df['LopNr'].apply(str)
    df[ei_col] = df[ei_col].apply(common.str2list)
    df[ci_col] = df[ci_col].apply(common.str2list)
    return df

def readgt():
    name = gs.names.treat_file
    return read_final(name)


def readtotal():
    list_cols = [
        'event_time_death',
	'event_time_inc',
	'event_time_nic',
	'surgery',
	'drug_class',
	'radio'
    ]

    name = gs.names.total_file
    df = read_final(name)
    f = lambda x: [] if x == 0 else x
    for c in list_cols:
        df[c] = df[c].apply(str).apply(common.str2list).apply(f)
    return df

"""
def readtimediff():
    name = gs.names.timediff_file
    return read_final(name)
"""


def readtimediff():
    # no type conversions but for LopNr
    name = gs.names.timediff_file
    df = read_generated_file(name)
    df['LopNr'] = df['LopNr'].apply(str)
    return df

