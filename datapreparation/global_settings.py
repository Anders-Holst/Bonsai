#! /home/jan/anaconda3/bin/python



""" -------------------------------

    Copyright (C) 2019 RISE
    This code was produced by RISE
    The 2019-08-29 version

    bonsai/src/global_settings.py

------------------------------------"""



import places


class long_names:
    input_data_group = 'group'
    dia_group     = 'Diagnostic_Group'
    out_p1        = 'outcome_p1_2019_12_05'
    out_p2        = 'outcome_p2_2019_10_16'
    outcome_file  = 'outcome_2019_12_05'
    timediff_file = 'timediff_2019_12_09'
    treat_file    = 'treat_2019_05_17'
    total_file    = 'total_2020_05_06'

class short_names:
    input_data_group = 'group'
    dia_group = 'group'
    outcome_file = 'outcome'
    timediff_file = 'timediff'
    treat_file = 'treat'
    total_file = 'total'


names = long_names



outcome_timegap    = 0 # 1826   # (nr of days in 5 years)
mean_diagnosis_age = 9.42 # (years, converted to days: d = y * 365.25)


diaspan  = ['1970-02-15', '2015-12-29']
incspan  = ['1970-01-12', '2016-12-30'] # sluten vard, changed 2019-12-04 
lifespan = ['1952-01-02', '2017-12-15']
nicspan  = ['1997-01-02', '2016-12-31']
drugspan = ['2005-06-28', '2017-12-31']

future   = '5000-01-01'   # any future date


"""
incare span:  1970-01-12 	 2016-12-30
nicare span:  1997-01-02 	 2016-12-31
dia span:     1970-02-15 	 2015-12-29
oc span:      1961-12-02 	 2024-12-01
life span:    1975-11-02 	 2017-12-15
drug span:    2005-06-28 	 2017-12-31
"""
