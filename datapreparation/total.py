#! /home/jan/anaconda3/bin/python



""" -------------------------------

    Copyright (C) 2018 RISE
    This code was produced by RISE
    The 2019-08-29 version

    bonsai/src/total.py


------------------------------------"""

# global imports

import pandas as pd
import numpy as np
from datetime import datetime as timeclass
from dateutil import parser
import copy

# local imports

import global_settings as gs
import bonsai_io as bio
import common
import merge
import diagnose
import bonsai_time as btime


# import time
import analyze


""" ----------------------------------------

    1 dia

---------------------------------------"""

def dia_base():
    dia = bio.read_generated_dia()
    grp = bio.read_diagroups()
    dia = diagnose.add_group(dia, grp)
    # dia = common.rmcols(dia, ['ICD10', 'SNOMED']) # real
    dia = common.rmcols(dia, ['SNOMED']) # real
    dia = diagnose.add_pers(dia)
    dia['DIA'] = 1
    return dia


""" ----------------------------------------

    2 control group and base set

---------------------------------------"""



def control_base():
    ctr = bio.readcontrol()
    ctr['DIA'] = 0
    m = gs.mean_diagnosis_age
    f = lambda x: common.add_years_to_full_year(x, m)
    ctr['DiagnosDat'] = ctr['FoddAr'].apply(f)
    return ctr


def base_set():
    dia = dia_base()
    ctr = control_base()

    b = gs.places.base_columns # real
    c = ctr.columns
    d = dia.columns

    """
    print('base 1')
    for z in c: print ('\t', z)
    print()
    print('base 2')
    for z in d: print ('\t', z)
    print()
    print('base 3')
    for z in b: print ('\t', z)
    print()
    """

    
    for col in list(set(b) - set(c)):
        ctr[col] = np.nan
    for col in list(set(b) - set(d)):
        dia[col] = np.nan

    
    

    ctr = ctr[b]
    dia = dia[b]

    return pd.concat([dia, ctr])



""" ----------------------------------------

    3 cause  

---------------------------------------"""

def add_cause(df):                     # 0 for empty list is ok!
    cols = gs.places.cause_added_cols  # 'ICD', 'life_events', 'DODSDAT'
    df_icd10 = df[['LopNr', 'ICD10']]  # The LopNr and the recidive code
    cause = bio.readcause(df_icd10)    # don't read the recidive class causes
    df = diagnose.add_cols(df, cause, 'LopNr', cols)

    c0 = 'DiagnosDat'
    c1 = 'DODSDAT'
    # df['Life(s)'] = df.apply(lambda d: date_diff_in_seconds(d, c1, c0), axis = 1)
    # df['Life(days)'] = df['Life(s)'] / (60*60*24)
    df['life'] = df.apply(lambda d: date_diff_in_days(d, c1, c0), axis = 1)

    return df


""" ----------------------------------------

    4 care 

---------------------------------------"""


def first_list(df1, df2, fcol, lex8, lex9, xcol = 'LopNr', ycol = 'DIAGNOS', zcol = 'INDATUM', split = True, n = 1,):
    df2 = btime.first_time_aho(df2, df1, xcol, ycol, zcol, fcol, lex8, lex9, split = split, n = n)
    return df2


def all_list(df1, df2, acol, lex8, lex9, xcol = 'LopNr', ycol = 'DIAGNOS', zcol = 'INDATUM', split = True, n = 1,):
    df2 = btime.all_times_aho(df2, df1, xcol, ycol, zcol, acol, lex8, lex9, split = split, n = n)
    return df2
  

def add_care(df, care, fcol, acol, n = 1, n2 = 0, n1 = 0):
    xcol = 'LopNr'
    
    if n2 > 0:
        care = care[n1:n2]                     # an initial part of incare

    lex8 = bio.read_care_icd8_lex()
    lex9 = bio.read_care_icd9_lex()

    
    #print('add ICD10 & DiagnosDat')
    #care = diagnose.add_cols_aho(care, df, xcol, ['ICD10','DiagnosDat'], True)

    L = list(df[xcol].values)
    care = care[ care[xcol].isin(L) ]      # part of care with LopNr in df

    print('    creating first_list')
    f = first_list(df, care, fcol, lex8, lex9, n = n)
        
    print('    adding first_list')
    df = common.add_list_cols(df, f, xcol, [fcol])
    
    print('    creating all_list')
    a = all_list(df, care, acol, lex8, lex9, n = n)
    
    print('    adding all_list')
    df = common.add_list_cols(df, a, xcol, [acol])
    return df   

def add_incare(df, n = 1, n2 = 0, n1 = 0):
    fcol = 'first_incare'
    acol = 'all_incare'
    inc  = bio.readincare()
    df   = add_care(df, inc, fcol, acol, n = n, n2 = n2, n1 = n1)
    return df


def add_nicare(df, n = 1, n2 = 0, n1 = 0):
    fcol = 'first_nicare'
    acol = 'all_nicare'
    nic  =  bio.readnicare()
    df   = add_care(df, nic, fcol, acol, n = n, n2 = n2, n1 = n1)
    return df

""" ----------------------------------------

    5 drugs

---------------------------------------"""


def add_drugs(df, n = 1, n2 = 0, n1 = 0): 
    xcol = 'LopNr'
    fcol = 'first_drug'
    acol = 'all_drug'
    ycol = 'ATC'
    zcol = 'EDATUM'
    drug  =  bio.readdrug(n)  # read n characters of the ATC code !!!

    print('drug.shape', drug.shape)
    if n2 > 0:
        drug = drug[n1:n2]                     # an initial part of drugs

    #ilex = [] 
    #if n == 1:
    #    ilex = bio.read_incare_lex() 
    

    # instead fo letting both first_list and all_list adding the columns
    # 'ICD10'  and 'DiagnosDat' it should be done here with
    # df = df1 and drug = df2
    
    L = list(df[xcol].values)
    drug = drug[ drug[xcol].isin(L) ]      # part of drug with LopNr in df

    #print('add ICD10 & DiagnosDat')
    #drug = diagnose.add_cols_aho(drug, df, xcol, ['ICD10','DiagnosDat'], True)

    print('    creating first_list')
    f = first_list(df, drug, fcol, False, False, xcol=xcol ,ycol=ycol, zcol=zcol, n = n)
    
    print('    creating all_list')
    a = all_list(df, drug, acol, False, False, xcol=xcol , ycol=ycol, zcol=zcol, n = n)
    
    print('    adding first_list')
    df = common.add_list_cols(df, f, xcol, [fcol])
    
    print('    adding all_list')
    df = common.add_list_cols(df, a, xcol, [acol])
        
    # drug.columns = ['LopNr', 'DIAGNOS', 'INDATUM']
    # drug = first_list(df, drug, fcol, split = False)
    # drug = first_list(df, drug, fcol, xcol=xcol , ycol=ycol, zcol=zcol)
    # df = common.add_list_cols(df, drug, xcol, [fcol])
    
    return df


""" ----------------------------------------

    care and cause merged

    not used in the present version

---------------------------------------"""



def concat(d, cs):
    L = []
    for c in cs:
        x = d[c]
        if common.isarray(x):
            L += list(x)
    return list(np.unique(L))


def merge_care_and_cause(df, diacol):
    nr = 0
    sv = 'DIA_SV'
    ov = 'DIA_OV'
    cc = 'life_events'
    cols = [sv, ov, cc]
    df = diagnose.add_incare(df, nr)
    df.columns = [sv if x=='DIAGNOS' else x for x in df.columns]
    df = diagnose.add_nicare(df, nr)
    df.columns = [ov if x=='DIAGNOS' else x for x in df.columns]
    df = add_cause(df)
    # df[diacol] = df.apply(lambda d:concat(d, cols), axis = 1)
    # df = common.rmcols(df, cols)
    # print(df.columns)
    return df

""" ----------------------------------------

    6 outcome

    The outcomes relates the input data the output data.
    The ouput is given by care, cause and drug
    outcome_01 give us the time of cause, if occurred and
    first time of inst. care for each class of cause and 
    cause diagnosis, as given by the ICD10 first letter.
    

    obs!!! the times produced by listdiff_*** are obsolet
    and not used and the assigments of base['Inc(s)'] and
    base['Inc(days)'] in outcome_01 may be removed.
    They are replaces by values given to df[ecol] and
    df[ne] in add_time_differences

---------------------------------------"""

def date_diff_in_seconds(df, c1, c0):
    return common.time_diff_in_seconds(df[c1], df[c0])

def date_diff_in_days(df, c1, c0):
    return common.time_diff_in_days(df[c1], df[c0])

def date_diff_in_years(df, c1, c0):
    return common.time_diff_in_years(df[c1], df[c0])

def list_diff_in_seconds(df, c1, c0):
    inc_list = df[c1]
    if not isinstance(inc_list, list):
        return inc_list
    L = []
    dtime = df[c0]
    for z in inc_list:
        if not isinstance(z, list):
            print (z)
        else:
            diff = common.time_diff_in_seconds(z[1], dtime)
            L = L + [[z[0], diff]]
    return L

def list_diff_in_days(df, c1, c0):
    inc_list = df[c1]
    if not isinstance(inc_list, list):
        return inc_list
    L = []
    dtime = df[c0]
    for z in inc_list:
        if not isinstance(z, list):
            print (z)
        else:
            diff = common.time_diff_in_days(z[1], dtime)
            L = L + [[z[0], diff]]
    return L



def outcome(n = 1, save = False, n2 = 0, n1 = 0):
    # n is the number of considered characters of the icd10 code 

    name = gs.names.outcome_file
    base = base_set()
    base = add_cause(base)  # recidiv not added    
    print('\n    add incare\n') 
    base = add_incare(base, n = n, n2 = n2, n1 = n1) # recidiv not added
    print('\n    add nicare\n') 
    base = add_nicare(base, n = n, n2 = n2, n1 = n1) # recidiv not added    
    base = add_drugs(base, n = n, n2 = n2, n1 = n1)  # just a single letter for each order

    if save:
        bio.save_generated_file(base, name)

    return base

# instead of running all of the outcome we may divide it into parts
# storing and reading each partial computation


def out_p1(n =3, n2 = 0, n1 = 0):

    out_p1 = gs.names.out_p1
    base = base_set()
    base = add_incare(base, n = n, n2 = n2, n1 = n1)
    base = add_nicare(base, n = n, n2 = n2, n1 = n1) 
    bio.save_generated_file(base, out_p1)


# the drugs are not used at present and we may skip them

def out_no_drugs(n =3, n2 = 0, n1 = 0):
    out_p1 = gs.names.out_p1
    out_file = gs.names.outcome_file
    df = bio.read_generated_file(out_p1)
    df['LopNr'] = df['LopNr'].apply(str)
    df = add_cause(df)
    bio.save_generated_file(df, out_file)
   
# with the drugs we have three parts    

def out_p2(n =3, n2 = 0, n1 = 0):

    out_p1 = gs.names.out_p1
    out_p2 = gs.names.out_p2    
    df = bio.read_generated_file(out_p1)
    df['LopNr'] = df['LopNr'].apply(str)
    df = add_drugs(df, n = n, n2 = n2, n1 = n1)
    bio.save_generated_file(df, out_p2)

def out_p3(n =3, n2 = 0, n1 = 0):
    out_p2 = gs.names.out_p2
    out_file = gs.names.outcome_file
    df = bio.read_generated_file(out_p2)
    df['LopNr'] = df['LopNr'].apply(str)
    df = add_cause(df)  
    bio.save_generated_file(df, out_file)
 
    
""" ----------------------------------------

    7 time differences
    
---------------------------------------"""



def time_differences(df, span, fcol):
    
    gap = gs.outcome_timegap
    dia = df['DiagnosDat']
    
    [r1, r2] = span    
    p1       = common.add_days(dia, gap)
    p2       = df['DODSDAT']                
    p3       = df['last_date']
    f_str    = df[fcol]
    emigrated = df['VitalStatus'] == 'emigrated'

    start = common.time_max(r1, p1)
    ne_end = r2
    if p2 != '0':
        ne_end = common.time_min(p2, r2)
    elif emigrated and common.isstr(p3):   # zzz
        ne_end = common.time_min(p3, r2)
        

    no_event_time = common.time_diff_in_days(ne_end, start)    
    no_event_time = max(0, no_event_time)
    event_times = []

    if f_str == '[]' or f_str == '0':
        return [event_times, no_event_time]

    f_list = common.str2list(f_str)
    
    for event in f_list:
        [icd10_class, event_date] = event
        class_time = common.time_diff_in_days(event_date, start)   # cannot be negative !
        if class_time < 0: print('event_time < 0')
        event_times += [[icd10_class, class_time]]

    return [event_times, no_event_time]



def no_event_time_death(df):

    gap = gs.outcome_timegap
    dia = df['DiagnosDat']
    [r1, r2] = gs.lifespan                  # register start and end
    p1       = common.add_days(dia, gap)
    p2       = df['DODSDAT']                # we need a DODSDAT str it must have been convered
    p3       = df['last_date']
    emigrated = df['VitalStatus'] == 'emigrated'

    if not common.isstr(p2):
        print('DODSDAT need to be converted to a str')

    start = common.time_max (r1 , p1 )      # with a complete register it's p1
    end = r2
    if p2 != '0':
        end = p2

    elif emigrated and common.isstr(p3):   # zzz
        end = common.time_min(p3, r2)
    time = common.time_diff_in_days(end, start)
    return max(0, time)
    


def add_no_event_time_death(df, scol):     # no_event_time_death
    df[scol]   = df.apply(no_event_time_death, axis = 1)
    return df

def add_time_differences(df, span, ecol, fcol): # zzz
    ne = 'no_' + ecol 
    df['tmp']  = df.apply(lambda d: time_differences(d, span, fcol), axis = 1)
    df[ecol]  = df['tmp'].apply(lambda x: x[0])
    df[ne] = df['tmp'].apply(lambda x: x[1])
    df = df.drop(['tmp'], axis=1)
    return df

def all_timedifferences(df, span, acol):

    if isinstance(df[acol], str):
        a  =  common.str2list(df[acol])
    else:
        a = []
        print("Bad string: ", df[acol])
        
    if a == []:
        return []
    
    dia  =  df['DiagnosDat']
    gap  =  gs.outcome_timegap
    gap_date  = common.add_days(dia, gap)
    [start, end] = span
    z    = common.time_max(gap_date, start)
    data = []
    
    for x in a:
        # y = x[0]
        dates = x[1]
        dates.sort()
        diffs = []

        for date in dates:
            diff = common.time_diff_in_days(date, z)
            diffs += [diff]
        data += [[x[0], diffs]]
        
    return data



def add_all_timedifferences(df, span, dcol, acol):
    df[dcol] = df.apply(lambda d: all_timedifferences(d, span, acol), axis = 1)
    return df


def censored_time(df, spanstart):
    """
    censored time = max (0, min (r1 , p2 , p3 ) − p1 )
    """
    
    gap = gs.outcome_timegap
    dia = df['DiagnosDat']
    
    r1  = spanstart    
    p1  = common.add_days(dia, gap)
    p2  = df['DODSDAT']                
    p3  = df['last_date']
    emigrated = df['VitalStatus'] == 'emigrated'

    # start = p1
    end = r1
    if p2 != '0':
        end = common.time_min(p2, r1)

    elif emigrated and common.isstr(p3):   # zzz   
        end = common.time_min(p3, r1)
    
    time = common.time_diff_in_days(end, p1)

    return max(0, time)




def add_censored_times(df):
    [incstart, incend] = gs.incspan
    [nicstart, nicend] = gs.nicspan
    [drugstart, drugend] = gs.drugspan
    df['censored_time_inc']  = df.apply(lambda d: censored_time(d, incstart), axis = 1)
    df['censored_time_nic']  = df.apply(lambda d: censored_time(d, nicstart), axis = 1)
    df['censored_time_drug']  = df.apply(lambda d: censored_time(d, drugstart), axis = 1)
    df['censored_time_death']  = 0
    return df
 

def timediff():
    name = gs.names.outcome_file

    ispan = gs.incspan
    nspan = gs.nicspan
    dspan = gs.drugspan

    ficol = 'first_incare'
    fncol = 'first_nicare'
    fdcol = 'first_drug'

    eicol = 'event_time_inc'
    encol = 'event_time_nic'
    edcol = 'event_time_drug'

    aicol = 'all_incare'
    ancol = 'all_nicare'
    adcol = 'all_drug'

    dicol = 'diff_incare'
    dncol = 'diff_nicare'
    ddcol = 'diff_drug'    

    scol = 'no_event_time_death'
    df   = bio.read_generated_file(name)

    # 'DODSDAT' is read as an int
    df['DODSDAT'] = df['DODSDAT'].apply(str)

    # The three emigrated persons without a last date are removed
    df = df[(df['VitalStatus'] != 'emigrated') | df['last_date'].notnull()]

    df   = add_no_event_time_death(df, scol)
    df   = add_time_differences(df, ispan, eicol, ficol)
    df   = add_time_differences(df, nspan, encol, fncol)

    df   = add_censored_times(df)
    # df   = add_time_differences(df, dspan, edcol, fdcol)
    
    # df   = add_all_timedifferences(df, ispan, dicol, aicol)
    # df   = add_all_timedifferences(df, nspan, dncol, ancol)
    # df   = add_all_timedifferences(df, dspan, ddcol, adcol)
    

    #df   = df[gs.places.out_selection_01]
    df   = add_other_str(df)
    bio.save_generated_file(df, gs.names.timediff_file)



""" ----------------------------------------

    8 counts and frequences
    (replaced by 9)

---------------------------------------"""


def capitals():
    c = []
    for x in range(65, 91):
        c += [chr(x)]
    return c

def event_types(events):
    t = []
    for event in events:
        t = t + [event[0]]
    return set(t)

def event_data(df, nr = 0):

    ev_col  = 'event_time_inc'
    nev_col = 'no_event_time_inc'
    caps = set(capitals())
    ev  = {}
    nev = {}
    for c in  caps:
        ev[c]  = []
        nev[c] = []
    data = df[[ev_col, nev_col]].values 
    if nr > 0:
        data = data[:nr]
    for line in data:
        events   = line[0]
        nev_time = line[1]
        nev_types = caps -  event_types(events)
        for c in nev_types:
            nev[c] += [nev_time]
        for event in events:
            [c, ev_time] = event
            if c.isupper():
                ev[c] += [ev_time]
    return ev, nev


def extract(d):
    xd = {}
    for k in d:
        t = d[k]
        xd[k] = {}
        xd[k]['nr']   = len(t)
        xd[k]['days'] = sum(t)
    return xd

def xdata(df, nr = 0):
    ev, nev = event_data(df, nr = nr)
    xev  = extract(ev)
    xnev = extract(nev)

    total = {}
    for c in capitals():
        total[c] = {}
        nr       = xev[c]['nr'] 
        days     = xev[c]['days'] + xnev[c]['days'] 

        total[c]['nr']   = nr
        total[c]['days'] = days
        # total[c]['freq'] = round((nr / days) * int(1e6), 2)

    return total


def print_xd(xd):
    # print('type\t nr\t 1000 years\t nr / (1000 years)')
    print('type\t nr\t nr / (1000 years)')
    print()
    for c in capitals():
        nr = xd[c]['nr']
        y  = xd[c]['days'] / (365.24 * 1000)
        if y > 1e-10:
            f  = nr / y
        else:
            f = -1
        if nr > 0:
            print(c + '\t', nr, '\t', round(y, 1), '  \t', round(f,3))




""" ----------------------------------------

    9 frequencies and probabilities

    compare diagnosis group and control group by 
    the probability that the insitutional care 
    frequency being higher for the fist category
    
---------------------------------------"""


def compute_prob(xd, xc):

    sf_lim  = np.log(1e-7)
    data = []
    for c in capitals():
        nd = xd[c]['nr']
        nc = xc[c]['nr']
        sd = xd[c]['days']
        sc = xc[c]['days']
        if nd > 0 and nc > 0:
            fd = round(1e6 * nd/sd, 3)
            fc = round(1e6 * nc/sc, 3)
            sf  = common.exp2sample_logsf(nd, nc, sd, sc)
            H1  = int(sf < np.log(1e-7))
            nsf  = round(-sf, 2) 
            data = data + [ [c, nd, nc, sd, sc, fd, fc, nsf, H1] ]
    
    df = pd.DataFrame(data)
    df.columns = [
        'group',
        'nr_dia',
        'nr_ctr',
        'days_dia',
        'days_ctr',
        'freq_dia',
        'freq_ctr',
        '-log(1-p)',
        '1-p < 1e-7',
    ]
    return df


selected_cols = [
    'group',
    'nr_dia',
    'nr_ctr',
    'freq_dia',
    'freq_ctr',
    '-log(1-p)',
    '1-p < 1e-7',
    ]


def compare_inc_group(group = 0):
    group_col = 'Diagnostic_Group'
    td = bio.readtimediff()

    # td = td[td['oth'] == 0]  # remove all other diseases
    dia = td[td['DIA'] == 1]
    if group > 0:
        dia = dia[dia[group_col] == group]
    ctr = td[td['DIA'] == 0]
    xd = xdata(dia)
    xc = xdata(ctr)
    df = compute_prob(xd, xc) # lägg till styrka ?
    # df = df[selected_cols]
    df['M days_dia'] = df['days_dia']/1e6
    df['M days_ctr'] = df['days_ctr']/1e6
    df = common.rmcols(df, ['days_dia', 'days_ctr'])
    return df



def compare_inc():
    df =  compare_inc_group()
    bio.save_local_result(df, 'sluten_vard_5y_alla')
    for group in range(1, 13):
        df = compare_inc_group(group)
        bio.save_local_result(df, 'sluten_vard_5y_grupp_' + str(group))




""" ----------------------------------------

    10 groups of other diseases

---------------------------------------"""



def selection(df, Take, Leave = []):
    annat = 'AnnatSpecDiagnos'
    a = df[annat]
    if not common.isstr(a):
        return np.nan
    for x in Leave:
        if (x in a):
            return False
    for x in Take:
        if (x in a):
            return True
    return False
      

 
def add_other(df):

    bw       = lambda d: selection(d, Take = ['Beck'])
    ep       = lambda d: selection(d, Take = ['pilep'])
    men      = lambda d: selection(d, Take = ['Sip', 'MEN'])
    men2     = lambda d: selection(d, Take = ['Sip', 'MEN'], Leave = ['1'])
    nf       = lambda d: selection(d, Take =  ['eckling', 'NF', 'euro'])
    nf1      = lambda d: selection(d, 
                                    Take =  ['eckling', 'NF', 'euro'],
                                    Leave = ['typ 2', 'typ II', 'NF2', 'NF 2'])
    t21      = lambda d: selection(d, Take = ['21', 'own'])
    tub      = lambda d: selection(d, Take = ['Tub'])

    df['bw']   = df.apply(bw,  axis = 1)
    df['ep']   = df.apply(ep,  axis = 1)
    df['men']  = df.apply(men, axis = 1)
    df['nf']   = df.apply(nf1, axis = 1)
    df['t21']  = df.apply(t21, axis = 1)
    df['tub']  = df.apply(tub, axis = 1)
    df['spe']  = df['bw'] + df['ep'] + df['men'] + df['nf'] + df['t21'] 
    df['oth']  = ~df['AnnatSpecDiagnos'].isnull()

    return df 

def add_other_str(df):
    cols  = ['bw', 'ep', 'men', 'nf', 't21', 'tub', 'oth']
    df    = add_other(df)
    df    = common.df_bool2str(df, cols)
    return df 
    


""" ----------------------------------------

    11 surgery

---------------------------------------"""

def surg_codes(df):
    kod = df['Kod']
    sub = df['Sub']

    if common.isnull(kod):
        return np.nan

    kod = kod.replace('+', ' ')
    kod = kod.replace('/', ' ')
    kod = kod.split()

    if not common.isnull(sub):   

        sub = sub.replace('+', ' ')
        sub = sub.split()
        for k in kod:
            if not k[0] == sub[0][0]:
                sub += k
        return sub
    return kod    


def generate_surgery():
    df = bio.readsurgery()
    df = df.sort_values('LopNr').reset_index(drop=True)
    df['surgery'] = df.apply(surg_codes, axis = 1)
    df = common.compression_with_dates(df, 'LopNr', 'surgery', 'Datum')
    bio.save_generated_file(df, 'surgery')



def surgery_diff(df):
    dia     = df['DiagnosDat']
    surgery = df['surgery']
    surgeries_and_diffs = []
    for [codes, date] in surgery:
        diff = common.time_diff_in_days(date, dia)
        surgeries_and_diffs += [[codes, diff]]
    return common.code_compress(surgeries_and_diffs)

 
def add_surgery(df):
    surg = bio.read_generated_surgery()
    cols = list(surg.columns) 
    cols.remove('LopNr')
    df = common.add_list_cols(df, surg, 'LopNr', cols)
    df['surgery_diff']  = df.apply(surgery_diff, axis = 1)  
    return df



def surgery_ohe(df):
    c = 'kirurgi'
    df = diagnose.one_x(df, c)
    df = common.rmcols(df, [c])
    return df



""" ----------------------------------------

    12 cytostatica

---------------------------------------"""



def drug_classes(x, treatment, cyto_code, drug_class):

    tx = treatment[ treatment['LopNr'] == x ]
    if tx.empty:
        return []
    
    tx = tx[['BehNr', 'cyto_date']]
    tx = tx.sort_values('cyto_date')
    bs = np.unique(tx['BehNr'].values)  # all the 'BehNr' for x, max 2
    category_date_list = []
    
    for b in bs:
        codes = common.look_up_all_values(b,  cyto_code,  'BehNr', 'DrogKod')
        categories = []
        for c in np.unique(codes):
            categories += common.look_up_all_values(c,  drug_class,  'Drogkod', 'Kategori')
        categories = np.unique(categories)
        txb = tx[tx['BehNr'] == b]
        dates = np.unique(txb['cyto_date'].values)
        d0 = dates[0]
        for c in categories:
            category_date_list += [[c, d0]]
        """
        for d in dates:
            for c in categories:
                category_date_list += [[c, d]]
        """
        
    return category_date_list



def drug_diff(df):
    dia  = df['DiagnosDat']
    drugs = df['cytoclass']
    drugs_and_diffs = []
    for [code, date] in drugs:
        diff = common.time_diff_in_days(date, dia)
        drugs_and_diffs += [[code, diff]]
    return drugs_and_diffs
    
  

def add_cytostatica(df, treat):

    treat = treat[['LopNr', 'BehNr', 'cyto_date']]
    treat = treat.dropna()

    xs = np.unique(treat['LopNr'].values)
    print('    Nr of LopNr:s in treatment data with cyto', len(xs))
 
    
    # treat   = bio.readtreatment()   # LopNr > BehNr
    codes   = bio.readcytocodes()     # BehNr > DrogKod
    classes = bio.readdrugclasses()   # Drogkod > Kategori

    f = lambda x: drug_classes(x, treat, codes, classes)

    df['cytoclass'] = df['LopNr'].apply(f)
    df['cytoclass_diff']  = df.apply(drug_diff, axis = 1)  

    return df




""" ----------------------------------------

    13 stemmcell

---------------------------------------"""

"""
def stemcell(d, treat):
    x = d['LopNr']
    allo = common.look_up_all_values(x, treat, 'LopNr', 'Allogen')
    d['allo'] = 0
    d['auto'] = 0
    if 'Ja' in allo:
        d['allo'] = 1
    if 'Nej' in allo:
        d['auto'] = 1
    return d
"""

def stemcell(x, treat):
    tx = treat[ treat['LopNr'] == x ]
    if tx.empty:
        return [] 
    tx = tx[['Allogen', 'stem_date']]
    allo = tx[tx['Allogen'] == 'Ja']
    auto = tx[tx['Allogen'] == 'Nej']
    allo_dates = np.unique(allo['stem_date'].values)
    auto_dates = np.unique(auto['stem_date'].values)

    stem = []
    for date in allo_dates:
        stem += [['allo', date]]
    for date in auto_dates:
        stem += [['auto', date]]

    return stem


def stem_diff(df):
    dia  = df['DiagnosDat']
    stem = df['stemcell']
    stem_and_diffs = []
    for [code, date] in stem:
        diff = common.time_diff_in_days(date, dia)
        stem_and_diffs += [[code, diff]]
    return stem_and_diffs
    


def add_stemcell(df, treat):
    
    treat = treat[['LopNr', 'Allogen', 'stem_date']]
    treat = treat.dropna()
    xs = np.unique(treat['LopNr'].values)
    print('    Nr of LopNr:s in treatment data with stem', len(xs))
    df['stemcell'] = df['LopNr'].apply(lambda x: stemcell(x, treat))
    df['stemcell_diff']  = df.apply(stem_diff, axis = 1)  
    return df



def analyze_stemcell_data(df):   # not used for production of total
    cols = gs.places.stemcell_columns 
    cell  = bio.readtreatment()
    cell  = cell[['LopNr'] + cols]
    cell  = cell.drop_duplicates()
    cell  = cell[ ~( cell['Allogen'].isnull() & cell['Autolog'].isnull()) ]
    print(cell.shape)    
    print(cell[0:20]) 

    analyze.compare_lopnr(df, cell, 'base', 'cell')
    
    df   = diagnose.add_cols(df, cell, 'LopNr', cols)
    return df


""" ----------------------------------------

    14 stralung

---------------------------------------"""

# def look_up_

def radio_datum(treat, x, b):
    pass

def radio_classes(x, treat, codes, classes):

    cx = codes[codes['LopNr'] == x]
    
    if cx.empty:
        return []
    
    cx = cx[['BehNr','Target_Kod']]
    tx = treat[treat['LopNr'] == x]
    bs = np.unique(cx['BehNr'].values)  # all the 'BehNr' for x, max 2

    category_date_list = []
    for b in bs:
        cxb = cx[cx['BehNr'] == b]
        txb = tx[tx['BehNr'] == b]
        dates = np.unique(txb['radio_date'].values)
        codes = np.unique(cxb['Target_Kod'].values)

        d0 = dates[0]

        for c in codes:
            for k in classes.columns[1:]:
                categories = common.look_up_nonnan(c, classes,  'Target_Kod', k)
                if not categories == []:
                    t = categories[0]
                    if [t, d0] not in category_date_list:
                        category_date_list += [[t, d0]]
        
        """
        for d in dates:  
            for c in codes:
                for k in classes.columns[1:]:
                    categories = common.look_up_nonnan(c, classes,  'Target_Kod', k)
                    if not categories == []:
                        t = categories[0]
                        if [t, d] not in category_date_list:
                            category_date_list += [[t, d]]
        """

    return category_date_list


def radio_diff(df):
    dia  = df['DiagnosDat']
    radio = df['stralung']
    radio_and_diffs = []
    for [code, date] in radio:
        diff = common.time_diff_in_days(date, dia)
        radio_and_diffs += [[code, diff]]
    return radio_and_diffs
    


def add_radio(df, treat):

    treat = treat[['LopNr', 'BehNr', 'radio_date']]
    treat = treat.dropna()
    xs = np.unique(treat['LopNr'].values)
    print('    Nr of LopNr:s in treatment data with radio', len(xs))

    
    codes  = bio.readradiocodes()     # LopNr > Target_Kod (radioterapi)
    codes = codes.dropna()    
    classes = bio.readradioclasses()   # Target_kod > "Kategorier"
    f = lambda x: radio_classes(x, treat, codes, classes)
    df['stralung'] = df['LopNr'].apply(f)
    
    df['radio_diff']  = df.apply(radio_diff, axis = 1)  

    return df




""" ----------------------------------------

    14 total

---------------------------------------"""


def stem(df):
    z = []
    if df['allo'] == 1: z += ['allo']
    if df['auto'] == 1: z += ['auto']    
    return z

def other(df):
    z = []
    for s in ['bw', 'ep', 'men', 'nf', 't21', 'tub']:	
        if df[s] == '1': z += [s]
    return z


def event_time_death(df):
    events = df['life_events']
    if events == '0':
        return []

    events = common.str2list(events)

    gap   = gs.outcome_timegap
    dia   = df['DiagnosDat']
    start = common.add_days(dia, gap)
    end   = str(df['DODSDAT'])
    time  = common.time_diff_in_days(end, start)
    
    # t= int(df['life']) # redo this
    z = []
    for e in events:
        z += [[e, time]]
    
    return z
 
def kon(x):
    if x == '1': return 'pojke'
    if x == '2': return 'flicka'
    return x


def place_in_order(s, y, z):
    # place z after y in the list s 
    t = []
    for x in s:
        if not x == z:
            t += [x]
        if x == y:
            t += [z]
    return t



cols_to_remove = [
    'AnnatSpecDiagnos',
    'Diagnostic_Group',
    'K_n',
    'Diagnos_lder',
#    'DiagnosDat',
    'SlutBehDat',
    'FoddAr',
    'first_incare',
    'first_nicare',
    'all_incare',
    'all_nicare',
    'first_drug',
    'all_drug',
    'surgery',
    'cytoclass',
    'stemcell',
    'radio',
    'ICD',
    'DODSDAT',
    'life_events',
    'life',
    'oth',
    'bw',
    'ep',
    'men',
    'nf',
    't21',
    'tub',
    'spe',
#    'allo',
#    'auto',
]



def diagnosis_class(x):
    if isinstance(x, str):
        return x[:3]
    return x

def modify_total(df):
    
    # 1 change column names
    

    df = common.change_col_names(df, 'DIA', 'has_dia')
    # df = common.change_col_names(df, 'oth', 'has_other_dia') # other is removed
    df = common.change_col_names(df, 'stralung', 'radio')
    df = common.change_col_names(df, 'ICD10', 'diagnosis')
    

    # 2 define other_dia columns
     
    dfcopy = copy.copy(df)  # alt. = copy.deepcopy(df)
    # dfcopy['stemcell']  = df.apply(lambda d: stem(d), axis = 1)
    dfcopy['other_dia'] = df.apply(lambda d: other(d), axis = 1)    
    df = dfcopy

    # 3 define event_time_death

    df['event_time_death'] = df.apply(lambda d: event_time_death(d), axis = 1)
    df['all_event_time_death'] = df['event_time_death'].apply(lambda lst: [[x[0], [x[1]]] for x in lst])

    # 4 add diagnosis class

    df['diagnosis_class'] = df['diagnosis'].apply(diagnosis_class)

    # 5 in col. K_n replace 1 with pojke and replace 2 with flicka
    
    df['sex'] = df['K_n'].apply(kon)
    
    # 6 rearrange columns

    c = list(df.columns)
    c = place_in_order(c, 'Diagnostic_Group', 'sex')
    c = place_in_order(c, 'DODSDAT', 'event_time_death')
    c = place_in_order(c, 'diagnosis', 'diagnosis_class')
    # c = place_after('diff_nic', c, 'no_event_time_nic')
    df = df[c]

    
    # 7 remove not to be used columns
    
    df = common.rmcols(df, cols_to_remove)

    return df


def generate_total():
    df = bio.readtimediff()
    treat = bio.readtreatment()
    df = add_surgery(df)
    df = add_cytostatica(df, treat)
    df = add_stemcell(df, treat)
    df = add_radio(df, treat)
    df = modify_total(df)
    bio.save_generated_file(df, gs.names.total_file)
    
    print('total columns')
    for c in df.columns:
        print('\t', c)
    print()

    t = df[['DiagnosDat', 'cytoclass', 'cytoclass_diff', 'stemcell', 'stemcell_diff', 'radio', 'radio_diff']]
    bio.save_tmp_file(t, 'total_dsr_2020_05_06')

    t = df[['surgery_diff', 'cytoclass_diff', 'stem_diff', 'radio_diff']]
    bio.save_tmp_file(t, 'total_treatment_2020_05_06')
