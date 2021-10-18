import sys
import time
import pandas as pd
import numpy as np
from datetime import datetime as timeclass
from dateutil import parser
import copy

import global_settings as gs
import bonsai_io as bio
import diagnose
import separate_tasks as sep
import total
import analyze
import lexicon
import common
import merge
import diagnose
import bonsai_time as btime


dia = diagnose.generate_dia()
dia = diagnose.add_pers(dia)
dia['DIA'] = 1
dia.drop_duplicates('LopNr',keep='first',inplace=True)

ctr = total.control_base()
b = gs.places.base_columns
c = ctr.columns
d = dia.columns
for col in list(set(b) - set(c)):
        ctr[col] = np.nan

for col in list(set(b) - set(d)):
        dia[col] = np.nan

base = pd.concat([dia[b], ctr[b]])
base.set_index('LopNr',drop=False,inplace=True)


base = total.add_cause(base)

inc  = bio.readorigincare()
nic  = bio.readorignicare()

base = total.add_care(base, inc, 'first_incare', 'all_incare', n = 3, n2 = 0, n1 = 0)
base = total.add_care(base, nic, 'first_nicare', 'all_nicare', n = 3, n2 = 0, n1 = 0)
base = total.add_drugs(base, n=4)


base['DODSDAT'] = base['DODSDAT'].apply(str)
base = base[(base['VitalStatus'] != 'emigrated') | base['last_date'].notnull()]
base['K_n']['5251'] = 'pojke'

base   = total.add_no_event_time_death(base, 'no_event_time_death')

for col in ['first_incare','all_incare','first_nicare','all_nicare','first_drug','all_drug','life_events']:
  base[col] = base[col].apply(str)

def fixdate00(t):
  if t[-2:]=='00':
    t = t[:-2] + "15"
  return t

base['DODSDAT'] = base['DODSDAT'].apply(fixdate00)

base   = total.add_time_differences(base, gs.incspan, 'event_time_inc', 'first_incare')
base   = total.add_time_differences(base, gs.nicspan, 'event_time_nic', 'first_nicare')
base   = total.add_time_differences(base, gs.drugspan, 'event_time_drug', 'first_drug')

base = total.add_all_timedifferences(base, gs.incspan, 'all_event_time_inc', 'all_incare')
base = total.add_all_timedifferences(base, gs.nicspan, 'all_event_time_nic', 'all_nicare')
base = total.add_all_timedifferences(base, gs.drugspan, 'all_event_time_drug', 'all_drug')

base = total.add_censored_times(base)
base = total.add_other_str(base)

gs.places.surgery = gs.places.local_input + 'kodad_kirurgi_20200318.csv'
sur=bio.readsurgery()
sur = sur.sort_values('LopNr').reset_index(drop=True)
sur['surgery'] = sur.apply(total.surg_codes, axis = 1)
sur = common.compression_with_dates(sur, 'LopNr', 'surgery', 'Datum')
cols = list(sur.columns) 
cols.remove('LopNr')
base = common.add_list_cols(base, sur, 'LopNr', cols)
base['surgery_diff']  = base.apply(total.surgery_diff, axis = 1)

treat = bio.readtreatment()
base = total.add_cytostatica(base, treat)
base = total.add_stemcell(base, treat)
base = total.add_radio(base, treat)

base = total.modify_total(base)

bio.save_generated_file(base, "total_210522")
