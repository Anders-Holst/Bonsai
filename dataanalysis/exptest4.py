"""
Copyright (C) 2018-2021 RISE Research Institute of Sweden AB

File: exptest4.py

Author: anders.holst@ri.se

"""

import pandas as pd
from math import *
from scipy.special import hyp1f1

from hist import *
from visual import *

eps = 0.001

# Brute force summation, in logarithms to avoid overflow
def log_hyper_1f1_negintab(a, b, x):
    m = -a
    lr = 0.0
    lhg = 0.0
    for k in range(m):
        lr += log(x*(a+k)/((k+1)*(b+k)))
        lhg += log(1 + exp(lr - lhg))
    return lhg

# r *= (x*(a+k) / ((k+1)*(b+k)))
# hg += r

def log_hyper_1f1(a, b, x):
    # b ej negativt heltal -> log(hyp1f1)
    # b neg heltal -> interpolera b+-eps
    eps = 0.0001
    if b<0 and b==int(b):
        ret = (hyp1f1(a, b-eps, x) + hyp1f1(a, b+eps, x))/2
    else:
        ret = hyp1f1(a, b, x)
    if isinf(ret):
        print("overflow at hyp1f1(%f, %f, %f)" % (a,b,x))
        return -inf
    elif ret<=0:
        print("negative hyp1f1(%f, %f, %f)" % (a,b,x))
        return -inf
    return log(ret)

def log_hyper_1f1_interpol(a, b, x):
    # a ej negativt heltal -> interpolera log(hyp1f1)
    if a==int(a):
        return log_hyper_1f1_negintab(int(a), b, x)
    else:
        r1 = log_hyper_1f1_negintab(floor(a), b, x)
        r2 = log_hyper_1f1_negintab(ceil(a), b, x)
        prop = a - floor(a)
        return r2*prop + r1*(1.0-prop)

# pa = (n, t)
def LogLambdaProb_un(pa1, ll):
    (n1, t1) = pa1
    if ll <= 0:
        return -inf if n1 > 0 else 0.0
    return n1 * log(ll) - t1*ll

def LogLambdaProb_pr_un(pa1, ll, pr):
    (n1, t1) = pa1
    (n0, t0) = pr
    if ll <= 0:
        return -inf if n1+n0 > 0 else 0.0
    return (n1+n0) * log(ll) - (t1+t0)*ll

def LogDifflambdaProbOne_un(pa1, pa2, dl):
    (n1, t1) = pa1
    (n2, t2) = pa2
    if dl == 0:
        ret = 0.0
    elif dl > 0.0:
        ret = (-dl * t2) + log_hyper_1f1_negintab(-n2, -n1 - n2, dl*(t1 + t2))
    else:
        ret = (dl * t1) + log_hyper_1f1_negintab(-n1, -n1 - n2, -dl*(t1 + t2))
    return ret

def LogDifflambdaProbOne_pr_un_old(pa1, pa2, dl, pr):
    (n1, t1) = pa1
    (n2, t2) = pa2
    (n0, t0) = pr
    if dl == 0:
        ret = 0.0
    elif dl > 0.0:
        ret = (-dl*(t2+t0)) + log_hyper_1f1_interpol(-n2-n0, -n1-n2-2*n0, dl*(t1+t2+2*t0))
    else:
        ret = (dl*(t1+t0)) + log_hyper_1f1_interpol(-n1-n0, -n1-n2-2*n0, -dl*(t1+t2+2*t0))
    return ret

def LogDifflambdaProbOne_pr_un(pa1, pa2, dl, pr):
    (n1, t1) = pa1
    (n2, t2) = pa2
    if dl == 0 or t1+t2 == 0.0 or n1+n2 == 0:
        ret = 0.0
    else:
        if prior_style == 'partaverage': # separate priors per condition
            (n0, t0) = (pr[0], pr[0]*(t1+t2)/(n1+n2+pr[0]))
        elif prior_style == 'partaveragezero': # inget extra event
            (n0, t0) = (pr[0], pr[0]*(t1+t2)/(n1+n2))
        else:
            (n0, t0) = pr
        if dl > 0.0:
            ret = (-dl*(t2+t0)) + log_hyper_1f1_interpol(-n2-n0, -n1-n2-2*n0, dl*(t1+t2+2*t0))
        else:
            ret = (dl*(t1+t0)) + log_hyper_1f1_interpol(-n1-n0, -n1-n2-2*n0, -dl*(t1+t2+2*t0))
    return ret

def LogDifflambdaProbOne_num1_old(pa1, pa2, pr):
    (n1, t1) = pa1
    (n2, t2) = pa2
    (n0, t0) = pr
    # hitta gränser för pa1 och pa2 separat
    func1 = lambda ll: LogLambdaProb_pr_un(pa1, ll, pr)
    func2 = lambda ll: LogLambdaProb_pr_un(pa2, ll, pr)
    (rng, vals) = find_calc_hrange(func1, (n1+n0+1)/(t1+t0), (n1+n0+1)/(5*(t1+t0)), 0.1, 25)
    mx1 = max(vals)
    a1 = rng[0]
    b1 = rng[-1]
    d1 = (b1 - a1)/(len(rng)-1)
    (rng, vals) = find_calc_hrange(func2, (n2+n0+1)/(t2+t0), (n2+n0+1)/(5*(t2+t0)), 0.1, 25)
    mx2 = max(vals)
    a2 = rng[0]
    b2 = rng[-1]
    d2 = (b2 - a2)/(len(rng)-1)
    return (mx1, a1, b1, d1, func1, mx2, a2, b2, d2, func2)

def LogDifflambdaProbOne_num1(pa1, pa2):
    (n1, t1) = pa1
    (n2, t2) = pa2
    if t1+t2==0.0 or n1+n2==0:
        return ()
    # hitta gränser för pa1 och pa2 separat
    func1 = lambda ll: LogLambdaProb_un(pa1, ll)
    func2 = lambda ll: LogLambdaProb_un(pa2, ll)
    (rng, vals) = find_calc_hrange(func1, (n1+1)/t1, (n1+1)/(5*t1), 0.1, 25)
    mx1 = max(vals)
    a1 = rng[0]
    b1 = rng[-1]
    d1 = (b1 - a1)/(len(rng)-1)
    (rng, vals) = find_calc_hrange(func2, (n2+1)/t2, (n2+1)/(5*t2), 0.1, 25)
    mx2 = max(vals)
    a2 = rng[0]
    b2 = rng[-1]
    d2 = (b2 - a2)/(len(rng)-1)
    return (mx1, a1, b1, d1, func1, mx2, a2, b2, d2, func2)

def LogDifflambdaProbOne_num2(tup, dl):
    if not tup:
        return 0.0
    (mx1, a1, b1, d1, func1, mx2, a2, b2, d2, func2) = tup
    # sen gå igenom för fixt dl och summera
    dd = min(d1, d2)
    aa = max(a1, a2-dl)
    bb = min(b1, b2-dl)
    sm = 0.0
    lst = []
    x = aa
    while x<bb:
        lst.append(func1(x) + func2(x+dl))
        x += dd
    if len(lst)==0:
        y = func1((aa+bb)/2) + func2((aa+bb)/2 + dl)
        return y - mx1 - mx2
    mx = max(lst)
    sm = sum(map(lambda y: exp(y-mx), lst))
    return log(sm) + mx - mx1 - mx2 if sm>0.0 else -inf

#def LogDifflambdaProb_un(da1, da2, dl):
#    return sum(list(map(lambda pa1,pa2: LogDifflambdaProbOne_un(pa1, pa2, dl), da1, da2)))

#def LogDifflambdaProb_pr_un(da1, da2, dl, pr):
#    pr = (pr[0]/len(da1), pr[1]/len(da1))
#    return sum(list(map(lambda pa1,pa2: LogDifflambdaProbOne_pr_un(pa1, pa2, dl, pr), da1, da2)))

def LogDifflambdaProb_num1(da1, da2):
    return list(map(lambda pa1,pa2: LogDifflambdaProbOne_num1(pa1, pa2), da1, da2))

def LogDifflambdaProb_num2(tupl, dl):
    return sum(list(map(lambda tup: LogDifflambdaProbOne_num2(tup, dl), tupl)))

prior_style = 'partaveragezero'

def set_prstyle(st):
    global prior_style
    prior_style = st

def DifflambdaHist(da1, da2, pa1, pa2):
    (n1, t1) = pa1
    (n2, t2) = pa2
    if t1 == 0.0 or t2 == 0.0 or (n1 == 0 and n2 == 0):
        return False
    tupl = LogDifflambdaProb_num1(da1, da2)
    func = lambda dl: LogDifflambdaProb_num2(tupl, dl)
    (xres, yres) = find_calc_hrange(func, 0.0, (n1+n2+1)/(t1+t2), 0.1, 20)
    return [xres, normalize_logprob(yres)]

def LambdaHist(pa1):
    (n1, t1) = pa1
    if t1 == 0.0 or n1 == 0:
        return False
    func = lambda dl: LogLambdaProb_un(pa1, dl)
    (xres, yres) = find_calc_hrange(func, (n1+1)/t1, (n1+1)/(5*t1), 0.2, 20)
    return [xres, normalize_logprob(yres)]

def normalize_logprob(vec):
    mx = max(vec)
    sm = 0.0
    for i in range(len(vec)):
        sm += exp(vec[i] - mx)
    norm = log(sm) + mx
    return [exp(x-norm) for x in vec]

def find_calc_hrange(func, start, step, mindiff, maxdiff):
    def zpush(v):
        return 0.0 if abs(v)<2e-11 else v
    fa = nan
    fb = nan
    a = nan
    b = nan
    m1 = zpush(start - 0.5*step)
    fm1 = func(m1)
    m2 = zpush(start + 0.5*step)
    fm2 = func(m2)
    left = []
    right = []
    while True:
        #print(a, fa, m1, fm1, m2, fm2, b, fb)
        if fm1 > fm2:
            if not isnan(b):
                right = [(b, fb)] + right
            (m, fm) = (m1, fm1)
            (b, fb) = (m2, fm2)
        else:
            if not isnan(a):
                left = [(a, fa)] + left
            (a, fa) = (m1, fm1)
            (m, fm) = (m2, fm2)
        if isnan(a):
            m1 = zpush(3*m - 2*b)
            fm1 = func(m1)
            (m2, fm2) = (m, fm)
        elif isnan(b):
            (m1, fm1) = (m, fm)
            m2 = zpush(3*m - 2*a)
            fm2 = func(m2)
        else:
            tmp = [fa, fm, fb]
            mx = max(tmp)
            mn = min(tmp)
            if mx-mn < mindiff or (b-a)<step/1024:
                break
            if (m-a) > (b-m)*1.5:
                dir = 0
            elif (m-a)*1.5 < (b-m):
                dir = 1
            else:
                dir = 0 if fa > fb else 1
            if dir == 0:
                m1 = zpush(0.5*(a+m))
                fm1 = func(m1)
                (m2, fm2) = (m, fm)
            else:
                (m1, fm1) = (m, fm)
                m2 = zpush(0.5*(m+b))
                fm2 = func(m2)
    delta = min(m-a, b-m)
    left = [(a, fa)] + left
    right = [(b, fb)] + right
    if m-a < m-b:
        right = [(m, fm)] + right
    else:
        left = [(m, fm)] + left
    tr = mx - maxdiff
    #print("^", mn, mx, tr)
    resx = []
    resy = []
    leftind = 0
    if delta < 1e-10:
        delta = 1e-3
    while leftind < len(left) and left[leftind][1] > tr:
        resx = [left[leftind][0]] + resx
        resy = [left[leftind][1]] + resy
        a = zpush(left[leftind][0] - delta)
        leftind += 1
        while leftind >= len(left) or a > left[leftind][0] + delta/2:
            fa = func(a)
            #print("L", a, fa)
            if isnan(fa) or isinf(fa) or fa < tr:
                break
            resx = [a] + resx
            resy = [fa] + resy
            a = zpush(a-delta)
    rightind = 0
    while rightind < len(right) and right[rightind][1] > tr:
        resx = resx + [right[rightind][0]]
        resy = resy + [right[rightind][1]]
        b = zpush(right[rightind][0] + delta)
        rightind += 1
        while rightind >= len(right) or b < right[rightind][0] - delta/2:
            fb = func(b)
            #print("R", b, fb)
            if isnan(fb) or isinf(fb) or fb < tr:
                break
            resx = resx + [b]
            resy = resy + [fb]
            b = zpush(b+delta)
    #print(a, m, b)
    return (resx, resy)

def selfunc(entry, val):
    if type(entry) == list:
        if len(entry)>0 and type(entry[0])==list:
            entry = list(map(lambda x:x[0], entry))
        if type(val) == list:
            for x in val:
                if x in entry:
                    return True
            return False
        else:
            return val in entry
    else:
        if type(val) == list:
            return entry in val
        else:
            return val == entry

def dictincr(dic, key, incr=1):
    if key not in dic:
        dic[key] = incr
    else:
        dic[key] += incr

def selassoc(sel, key):
    if sel is False:
        return []
    else:
        return [x[-1] for x in sel if key == x[0] or (len(x)==3 and key == x[1])]

def assoc(lst, key, defl):
    for ele in lst:
        if ele[0]==key:
            return ele[1]
    return defl

def assoc_a0(lst, key, defl):
    if type(lst)==list:
        for ele in lst:
            if type(ele)==list:
                if ele[0]==key:
                    return ele[1]
            else:
                if ele==key:
                    return 0
        return defl
    else:
        return 0 if lst==key else defl        

def subsets(lst, ord):
    if ord==0:
        return [[]]
    elif  ord > len(lst):
        return []
    else:
        rest = subsets(lst[1:], ord-1)
        return list(map(lambda l:[lst[0]]+l, rest)) + subsets(lst[1:],ord)

def issubset(l1, l2):
    for e in l1:
        if not e in l2:
            return False
    return True

def tupleadd(t1, t2):
    return tuple((t1[i]+t2[i] for i in range(min(len(t1),len(t2)))))
    #return (t1[0]+t2[0],t1[1]+t2[1])

def tuplescale(t1, scale):
    return tuple((t1[i]*scale for i in range(min(len(t1),len(t2)))))
    #return (t1[0]*scale,t1[1]*scale)

def listminus(l1, l2):
    return [e for e in l1 if not e in l2]

def read_total_file(file):
    df = pd.read_csv(file, sep='\t')
    for c in df:
        if type(df[c][0])==str and df[c][0][0] == '[':
            df[c] = df[c].apply(eval)
    return df

def valsincolumn_old(df, col):
    d = df[col]
    res = []
    for val in d:
        if type(val)==list:
            for e in val:
                if not e in res:
                    res.append(e)
        else:
            if not val in res:
                res.append(val)
    return sorted(res)

def valsincolumn(df, col):
    d = df[col]
    res = []
    for val in d:
        if type(val)==list:
            for e in val:
                if type(e)==list:
                    if not e[0] in res:
                        res.append(e[0])
                else:
                    if not e in res:
                        res.append(e)
        else:
            if not val in res and not (type(val)==float and isnan(val)):
                res.append(val)
    return sorted(res)

#--------------------------

# Ny approach för att undvika betingad utspädning: För varje
# seneffekt, i lämplig submängd av data (patienter, kön, etc) kolla om
# en viss faktor (egenskaper, diagnos, behandling) är signifikant
# korrelerad. Sen betinga på varje annan faktor som också är
# korrelerad, och se om signifikansen försvinner helt, alt frekvensen
# förändras signifikant. Lista de som inte försvinner/ändras, i
# signifikansordning. Håll reda på om alla försvinner vid korsvis
# betingning, eller om flera finns kvar. Helst ska varje signifikant
# korrelation förklaras av en faktor, alt flera faktorer av samma typ.

def select_data(df, sel1, sel2, selg):
    # sel = [(col, val)]    selg = [('otherdia', 't21')] [('DIA',1)] [False, 'stralung','ANSIKTE']
    # cond = [col1, col2, ] 
    if df.empty:
        return (df,df)
    if selg:
        test = df[df.columns[0]].apply(lambda entry: True)
        for ele in selg:
            if len(ele) == 3:
                test = test & df[ele[1]].apply(lambda entry: not selfunc(entry, ele[2]))
            else:
                test = test & df[ele[0]].apply(lambda entry: selfunc(entry, ele[1]))
        df = df[test]
    test = df[df.columns[0]].apply(lambda entry: True)
    for ele in sel1:
        if len(ele) == 3:
            test = test & df[ele[1]].apply(lambda entry: not selfunc(entry, ele[2]))
        else:
            test = test & df[ele[0]].apply(lambda entry: selfunc(entry, ele[1]))
    d1 = df[test]
    if sel2 is False:
        d2 = df[~test]
    else:
        test = df[df.columns[0]].apply(lambda entry: True)
        for ele in sel2:
            if len(ele) == 3:
                test = test & df[ele[1]].apply(lambda entry: not selfunc(entry, ele[2]))
            else:
                test = test & df[ele[0]].apply(lambda entry: selfunc(entry, ele[1]))
        d2 = df[test]
    return (d1, d2)

def select_data_alt(df, selalt, selg):
    if df.empty:
        return (df,df)
    if selg:
        test = df[df.columns[0]].apply(lambda entry: True)
        for ele in selg:
            if len(ele) == 3:
                test = test & df[ele[1]].apply(lambda entry: not selfunc(entry, ele[2]))
            else:
                test = test & df[ele[0]].apply(lambda entry: selfunc(entry, ele[1]))
        df = df[test]
    test = df[df.columns[0]].apply(lambda entry: False)
    for ele in selalt:
        if len(ele) == 3:
            test = test | df[ele[1]].apply(lambda entry: not selfunc(entry, ele[2]))
        else:
            test = test | df[ele[0]].apply(lambda entry: selfunc(entry, ele[1]))
    df = df[test]
    return df

def count_effect(df, eff, ncol, ecol):
    n = 0
    t = 0
    for i in range(len(df)):
        row = df.iloc[i]
        x = assoc(row[ecol], eff, False)
        if x is not False:
            n += 1
            t += x
        else:
            t += row[ncol]
    return (n, t)

def count_effect_mtag(df, eff, tags):
    ecols = ['event_time_' + tag for tag in tags]
    ncols = ['no_event_time_' + tag for tag in tags]
    #ccols = ['censored_time_' + tag for tag in tags]
    res = (0, 0, 0)
    for i in range(len(df)):
        row = df.iloc[i]
        x = min([assoc(row[ecol], eff, inf) for ecol in ecols])
        if x is not inf:
            res = tupleadd(res, (1, x, 1))
        else:
            x = min([row[ncol] for ncol in ncols])
            res = tupleadd(res, (0, x, 1))
    return res

#def count_cond_effect(df, eff, ncol, ecol, cols, ignore):
#    dic = {}
#    for i in range(len(df)):
#        row = df.iloc[i]
#        x = assoc(row[ecol], eff, False)
#        if x is not False:
#            incr = (1, x)
#        else:
#            incr = (0, row[ncol])
#        increment_tuple_value(row, dic, cols, ignore, incr)
#    return dic

def make_inc_nic_diag_dict(df):
    dic = {}
    for i in range(len(df)):
        nr=df['LopNr'][i]
        if nr not in dic:
            dic[nr] = set()
        dia=df['DIAGNOS'][i]
        if type(dia)==str:
            dic[nr].update(set(dia.split(' ')))
    return dic

def count_inc_nic_diag(df, dic, pair, selg, tags, eff):
    (dic1, dic2, dick) = get_cond_stats_daydiff_lopnr(df, pair, [], selg, tags, eff)
    resdic1 = {}
    resdic2 = {}
    for nr in dic1['']:
        if nr in dic:
            for diag in dic[nr]:
                if eff in diag:
                    dictincr(resdic1, diag)
    for nr in dic2['']:
        if nr in dic:
            for diag in dic[nr]:
                if eff in diag:
                    dictincr(resdic2, diag)
    return (resdic1, resdic2)

def get_cond_stats_daydiff_old(df, pair, conds, selg, tags, eff):
    # Vi får kuta igenom hela data, och för varje sample sortera in värden i rätt dict och key
    # Regeln är att ett villkor är sant om eff-dag är senare än ev cond-dag
    (dd, dummy) = select_data(df, [], False, selg)
    kl = [""]
    for cond in conds:
        kl = [(k+"_"+str(cond[-1]), k+"_~"+str(cond[-1])) for k in kl]
        kl = [k[0] for k in kl] + [k[1] for k in kl]
    dick = { k : kl.index(k) for k in kl}
    dic1 = { k : (0, 0, 0) for k in kl}
    dic2 = { k : (0, 0, 0) for k in kl}
    invk = { kl.index(k) : k for k in kl}
    bits = [2**i for i in range(len(conds))]
    ecols = ['event_time_' + tag for tag in tags]
    ncols = ['no_event_time_' + tag for tag in tags]
    for i in range(len(dd)):
        row = dd.iloc[i]
        x = min([assoc(row[ecol], eff, inf) for ecol in ecols])
        # here, sort it into correct key
        # för varje cond och pair, hitta senaste dag före eff-dag (minus marginal)
        ylst = [y if y<x-30 else -inf for y in [assoc_a0(row[col], val, -inf) for col,val in ([pair]+conds if pair else conds)]]
        y = max(max(ylst), 0)
        k = invk[sum([b if v==-inf else 0 for (b,v) in zip(bits,ylst if not pair else ylst[1:])])]
        td = dic2 if pair and ylst[0]==-inf else dic1
        if x is not inf:
            td[k] = tupleadd(td[k], (1, x-y, 1))
        else:
            x = min([row[ncol] for ncol in ncols])
            td[k] = tupleadd(td[k], (0, max(x-y, 0), 1))
    return (dic1, dic2, dick)

global_margin = 1826

def set_margin(marg):
    global global_margin
    global_margin = marg

def get_cond_stats_daydiff(df, pair, conds, selg, tags, eff):
    global global_margin
    # Vi får kuta igenom hela data, och för varje sample sortera in värden i rätt dict och key
    # Regeln är att ett villkor är sant om eff-dag är senare än ev cond-dag
    marg = global_margin
    (dd, dummy) = select_data(df, [], False, selg)
    kl = [""]
    for cond in conds:
        kl = [(k+"_"+str(cond[-1]), k+"_~"+str(cond[-1])) for k in kl]
        kl = [k[0] for k in kl] + [k[1] for k in kl]
    dick = { k : kl.index(k) for k in kl}
    dic1 = { k : (0, 0, 0) for k in kl}
    dic2 = { k : (0, 0, 0) for k in kl}
    invk = { kl.index(k) : k for k in kl}
    bits = [2**i for i in range(len(conds))]
    ecols = ['all_event_time_' + tag for tag in tags]
    ncols = ['no_event_time_' + tag for tag in tags]
    ccols = ['censored_time_' + tag for tag in tags]
    for i in range(len(dd)):
        row = dd.iloc[i]
        # get first event across all tags
        x,c = min([(assoc(row[ecol], eff, [inf])[0], row[ccol]) for ecol,ccol in zip(ecols,ccols)],key=lambda p:p[0]+p[1]) 
        # för varje cond och pair, hitta senaste dag före eff-dag (minus marginal)
        condlst = [y if y<x+c-marg else -inf if y > x+c else "Block" for y in [assoc_a0(row[col], val, -inf) for col,val in conds]]
        yval = assoc_a0(row[pair[0]], pair[1], -inf) if pair else 0
        td = dic2 if pair and yval==-inf else dic1
        yval = max(yval, 0)
        if "Block" in condlst or yval >= x+c-marg:
            continue
        # sort it into correct key
        k = invk[sum([b if v==-inf else 0 for (b,v) in zip(bits,condlst)])]
        if x is not inf:
            ymax = max(condlst + [yval]) if pair else max(condlst + [0])
            td[k] = tupleadd(td[k], (1, x-max(ymax+marg-c, 0), 1))
        else:
            x,c = min([(row[ncol],row[ccol]) for ncol,ccol in zip(ncols,ccols)], key=lambda p:p[0])
            td[k] = tupleadd(td[k], (0, max(x-max(marg-c, 0), 0), 1))
    return (dic1, dic2, dick)

def get_cond_stats_daydiff_lopnr(df, pair, conds, selg, tags, eff):
    (dd, dummy) = select_data(df, [], False, selg)
    kl = [""]
    for cond in conds:
        kl = [(k+"_"+str(cond[-1]), k+"_~"+str(cond[-1])) for k in kl]
        kl = [k[0] for k in kl] + [k[1] for k in kl]
    dick = { k : kl.index(k) for k in kl}
    dic1 = { k : [] for k in kl}
    dic2 = { k : [] for k in kl}
    invk = { kl.index(k) : k for k in kl}
    bits = [2**i for i in range(len(conds))]
    ecols = ['event_time_' + tag for tag in tags]
    ncols = ['no_event_time_' + tag for tag in tags]
    for i in range(len(dd)):
        row = dd.iloc[i]
        x = min([assoc(row[ecol], eff, inf) for ecol in ecols])
        ylst = [y if y<x-30 else -inf for y in [assoc_a0(row[col], val, -inf) for col,val in ([pair]+conds if pair else conds)]]
        y = max(max(ylst), 0)
        k = invk[sum([b if v==-inf else 0 for (b,v) in zip(bits,ylst if not pair else ylst[1:])])]
        td = dic2 if pair and ylst[0]==-inf else dic1
        if x is not inf:
            td[k].append(str(row['LopNr']))
    return (dic1, dic2, dick)

def get_cond_stats_one(df, pair, conds, selg, tags, eff):
    (d1, d2) = select_data(df, [pair], False, selg)
    kl = [""]
    dl1 = [d1]
    dl2 = [d2]
    for cond in conds:
        dl1 = [select_data(d1, [cond], False, []) for d1 in dl1]
        dl2 = [select_data(d2, [cond], False, []) for d2 in dl2]
        kl = [(k+"_"+str(cond[-1]), k+"_~"+str(cond[-1])) for k in kl]
        dl1 = [d1[0] for d1 in dl1] + [d1[1] for d1 in dl1]
        dl2 = [d2[0] for d2 in dl2] + [d2[1] for d2 in dl2]
        kl = [k[0] for k in kl] + [k[1] for k in kl]
    dic1 = { k : count_effect_mtag(d, eff, tags) for k,d in zip(kl,dl1)}
    dic2 = { k : count_effect_mtag(d, eff, tags) for k,d in zip(kl,dl2)}
    dick = { k : kl.index(k) for k in kl}
    return (dic1, dic2, dick)

def get_cond_stats_comb(df, conds, selg, tags, eff):
    (dd, dummy) = select_data(df, [], False, selg)
    kl = [""]
    dl = [dd]
    for cond in conds:
        dl = [select_data(d, [cond], False, []) for d in dl]
        kl = [(k+"_"+str(cond[-1]), k+"_~"+str(cond[-1])) for k in kl]
        dl = [d[0] for d in dl] + [d[1] for d in dl]
        kl = [k[0] for k in kl] + [k[1] for k in kl]
    dic = { k : count_effect_mtag(d, eff, tags) for k,d in zip(kl,dl)}
    datadic = { k : d for k,d in zip(kl,dl)}
    #dick = { k : kl.index(k) for k in kl}
    return (dic, datadic)

def estimatefactor(dlst1, dlst2):
    lf = 0.0
    sn = 0.0
    for ((n1,t1),(n2,t2)) in zip(dlst1,dlst2):
        if n1 > 0 and n2 > 0:
            lf += min(n1,n2)*(log(n2/t2) - log(n1/t1))
            sn += min(n1,n2)
    return exp(lf/sn) if sn>0 else 1.0

def dictolistwithprior(cdic1, cdic2, cdic0, tscale = 1.0):
    n1,t1,s1 = (0,0,0)
    n2,t2,s2 = (0,0,0)
    for k in cdic0:
        (tmp1, tmp2) = (cdic1[k], cdic2[k])
        n1 += tmp1[0]
        t1 += tmp1[1]*tscale
        s1 += tmp1[2]
        n2 += tmp2[0]
        t2 += tmp2[1]*tscale
        s2 += tmp2[2]
    if prior_style == 'noninfo':  # Noninformative prior
        prn = 1.0/len(cdic0)
        lst1 = [ (cdic1[k][0] - prn, cdic1[k][1]*tscale) for k in cdic0 ]
        lst2 = [ (cdic2[k][0] - prn, cdic2[k][1]*tscale) for k in cdic0 ]
    elif prior_style == 'average':  # Average prior
        prn = 1.0/len(cdic0)
        prt = prn*(t1+t2)/(n1+n2+1)
        lst1 = [ (cdic1[k][0] + prn, cdic1[k][1]*tscale + prt) for k in cdic0 ]
        lst2 = [ (cdic2[k][0] + prn, cdic2[k][1]*tscale + prt) for k in cdic0 ]
    elif prior_style == 'partaverage': # separate priors per condition
        prn = 1.0/len(cdic0)
        lst1 = []
        lst2 = []
        for k in cdic0:
            (nn1, tt1, ss1) = cdic1[k]
            (nn2, tt2, ss2) = cdic2[k]
            prt = prn*(tt1+tt2)/(nn1+nn2+prn)
            lst1.append((nn1+prn, (tt1+prt)*tscale))
            lst2.append((nn2+prn, (tt2+prt)*tscale))
    elif prior_style == 'partaveragezero': # inget extra event
        prn = 1.0/len(cdic0)
        lst1 = []
        lst2 = []
        for k in cdic0:
            (nn1, tt1, ss1) = cdic1[k]
            (nn2, tt2, ss2) = cdic2[k]
            if nn1+nn2 > 0:
                prt = prn*(tt1+tt2)/(nn1+nn2)
                lst1.append((nn1+prn, (tt1+prt)*tscale))
                lst2.append((nn2+prn, (tt2+prt)*tscale))
            else:
                lst1.append((nn1, tt1*tscale))
                lst2.append((nn2, tt2*tscale))
    else:  # No prior
        lst1 = [ (cdic1[k][0], cdic1[k][1]*tscale) for k in cdic0 ]
        lst2 = [ (cdic2[k][0], cdic2[k][1]*tscale) for k in cdic0 ]
    return (lst1, lst2, (n1,t1,s1), (n2,t2,s2))

def analysediff(cdic1, cdic2, cdic0):
    (dlst1, dlst2, (n1, t1, s1), (n2, t2, s2)) = dictolistwithprior(cdic1, cdic2, cdic0, 1.0/36525)

    hist = DifflambdaHist(dlst1, dlst2, (n1, t1), (n2, t2))
    fact = estimatefactor(dlst1, dlst2)
    if hist is not False:
        mean, var = hist_mean_var(hist)
        sig = hist_significance(hist, 0.0, mean)
    else:
        mean = 0.0
        var= 0.0
        sig = 1.0
    nn1 = sum(list(map(lambda p:p[0], dlst1)))
    tt1 = sum(list(map(lambda p:p[1], dlst1)))
    nn2 = sum(list(map(lambda p:p[0], dlst2)))
    tt2 = sum(list(map(lambda p:p[1], dlst2)))
    hist1 = LambdaHist((nn1, tt1))
    if hist1 is not False:
        mean1, var1 = hist_mean_var(hist1)
    else:
        mean1 = 0.0
    hist2 = LambdaHist((nn2, tt2))
    if hist2 is not False:
        mean2, var2 = hist_mean_var(hist2)
    else:
        mean2 = 0.0
    return {'mean':mean, 'std':sqrt(var), 'sig':sig, 'fact': fact, 'mean1':mean1, 'n1':n1, 't1':t1, 's1':s1, 'mean2':mean2, 'n2':n2, 't2':t2, 's2':s2, 'hist':hist, 'hist1':hist1, 'hist2':hist2}

def analysediff01(cdic1, cdic2, cdic0):
    (dlst1, dlst2, (n1, t1, s1), (n2, t2, s2)) = dictolistwithprior(cdic1, cdic2, cdic0, 1.0/36525)

    hist = DifflambdaHist(dlst1, dlst2, (n1, t1), (n2, t2))
    fact = estimatefactor(dlst1, dlst2)
    if hist is not False:
        mean, var = hist_mean_var(hist)
        sig = hist_significance(hist, 0.0, mean)
    else:
        mean = 0.0
        var= 0.0
        sig = 1.0
    nn1 = sum(list(map(lambda p:p[0], dlst1)))
    tt1 = sum(list(map(lambda p:p[1], dlst1)))
    nn2 = sum(list(map(lambda p:p[0], dlst2)))
    tt2 = sum(list(map(lambda p:p[1], dlst2)))
    mean1 = nn1/tt1 if tt1 > 0 else 0.0
    mean2 = nn2/tt2 if tt2 > 0 else 0.0
    return {'mean':mean, 'std':sqrt(var), 'sig':sig, 'fact': fact, 'mean1':mean1, 'n1':n1, 't1':t1, 's1':s1, 'mean2':mean2, 'n2':n2, 't2':t2, 's2':s2, 'hist':hist, 'hist1':False, 'hist2':False}

def analysediff0(cdic1, cdic2, cdic0):
    #dlst1 = cdictolist(cdic1, cdic0, 1.0/36525)
    #dlst2 = cdictolist(cdic2, cdic0, 1.0/36525)
    #if prior_style == 'noninfo':
    #    pr = (-1, 0) # Noninformative prior
    #elif prior_style in ['average','partaverage','partaveragezero']:
    #    n1 = sum(list(map(lambda p:p[0], dlst1)))
    #    t1 = sum(list(map(lambda p:p[1], dlst1)))
    #    n2 = sum(list(map(lambda p:p[0], dlst2)))
    #    t2 = sum(list(map(lambda p:p[1], dlst2)))
    #    pr = (1, 1*(t1+t2)/(n1+n2+1)) # Average prior
    #else:
    #    pr = (0, 0) # No prior
    #hist = DifflambdaHist(dlst1, dlst2, pr)
    (dlst1, dlst2, (n1, t1, s1), (n2, t2, s2)) = dictolistwithprior(cdic1, cdic2, cdic0, 1.0/36525)
    hist = DifflambdaHist(dlst1, dlst2, (n1, t1), (n2, t2))
    if hist is not False:
        mean, var = hist_mean_var(hist)
        sig = hist_significance(hist, 0.0, mean)
    else:
        mean = 0.0
        var= 0.0
        sig = 1.0
    return {'mean':mean, 'std':sqrt(var), 'sig':sig}

def analyseeffects1(df, selg, coldic, tags, eff):
    resdic = {}
    for col in coldic:
        for val in coldic[col]:
            #print((col,val))
            #(dic1, dic2, kdic) = get_cond_stats_one(df, (col, val), [], selg, tags, eff)
            (dic1, dic2, kdic) = get_cond_stats_daydiff(df, (col, val), [], selg, tags, eff)
            resdic[(col,val)] = analysediff(dic2, dic1, kdic)
    return resdic

def analyseeffects_back(df, selg, pair, tags, efflst):
    resdic = {}
    for eff in efflst:
        col,val = pair
        (dic1, dic2, kdic) = get_cond_stats_one(df, pair, [], selg, tags, eff)
        resdic[eff] = analysediff(dic2, dic1, kdic)
    return resdic

def analyseeffects2(df, selg, resdic1, sig, tags, eff):
    resdic2 = {}
    lst = []
    for pair in resdic1:
        if resdic1[pair]['sig'] <= sig:
            lst.append(pair)
    lst.sort(key=lambda x: resdic1[x]['sig'])
    for (ind,pair1) in reversed(list(enumerate(lst))):
        mnsig = 0.0
        mntmp = resdic1[pair1]
        for pair2 in lst:
            if pair1 != pair2:
                #print(pair1,pair2)
                #(dic1, dic2, kdic) = get_cond_stats_one(df, pair1, [pair2], selg, tags, eff)
                (dic1, dic2, kdic) = get_cond_stats_daydiff(df, pair1, [pair2], selg, tags, eff)
                tmp = analysediff(dic2, dic1, kdic)
                if tmp['sig'] > mnsig:
                    mnsig = tmp['sig']
                    tmp['cond'] = pair2
                    mntmp = tmp
        resdic2[pair1] = mntmp
        if mnsig > sig:
            del lst[ind]
    return resdic2

def movetosaved(remaining, saved, removed, dic):
    # alla som inte tas bort av andra än removed flyttas till saved
    changed = False
    for (ind,pair) in reversed(list(enumerate(remaining))):
        lst = dic[pair]
        ok = True
        for conds,sig in lst:
            ok = False
            for cond in conds:
                if cond in removed:
                    ok = True
            if not ok:
                break
        if ok:
            changed = True
            saved.append(pair)
            del remaining[ind]
    return changed

def movetoremoved(remaining, saved, removed, dic):
    # alla som tas bort av någon i saved flyttas till removed
    changed = False
    for (ind,pair) in reversed(list(enumerate(remaining))):
        lst = dic[pair]
        ok = False
        for conds,sig in lst:
            ok = True
            for cond in conds:
                if cond not in saved:
                    ok = False
            if ok:
                break
        if ok:
            changed = True
            removed.append(pair)
            del remaining[ind]
    return changed

def allexcept(lst, ele):
    return [e for e in lst if e != ele]

def analyseeffects2new(df, selg, resdic1, sig, tags, eff):
    ciidic = {}
    lst = []
    # välj ut dem med signifiant indirekt effekt
    for pair in resdic1:
        if resdic1[pair]['sig'] <= sig:
            lst.append(pair)
            ciidic[pair] = []
    # i första passet, betinga var och en på var och en av de andra,
    # spara lista på vilka som gör den insignifikant
    # ta bort dem som
    # 1) blir insignifikant av någon (som den själv inte gör insignifikant)
    # 2) och inte behövs för att göra någon annan insignifikant
    # I loop: spara dem som inte tas bort av nåt, släng dem som tas bort av dem,
    #  iterera med resten dvs spara av resten dem som inte tas bort av kvarvarande
    # i nästa pass betinga på par (och senare trippler) av kvarvarande
    # ta bort enligt analog princip
    remaining = lst
    saved = []
    totsaved = []
    for order in [1,2,3]:
        remaining = saved + remaining
        conds = subsets(remaining, order)
        for pair in lst:
            for cond in conds:
                if not pair in cond:
                    #(dic1, dic2, kdic) = get_cond_stats_one(df, pair, cond, selg, tags, eff)
                    (dic1, dic2, kdic) = get_cond_stats_daydiff(df, pair, cond, selg, tags, eff)
                    tmp = analysediff0(dic2, dic1, kdic)
                    if tmp['sig'] > sig:
                        ciidic[pair].append((cond,tmp['sig']))
        remaining = lst.copy()
        saved = []
        removed = []
        changed = True
        while changed:
            # alla som inte tas bort av andra än removed flyttas till saved
            # alla som tas bort av någon i saved flyttas till removed
            changed = False
            if movetosaved(remaining, saved, removed, ciidic):
                changed = True
            if movetoremoved(remaining, saved, removed, ciidic):
                changed = True
    if remaining:
        alternatives = []
        origsaved = saved.copy()
        origremaining = remaining.copy()
        akeys = []
        for pair in remaining:
            for sp in ciidic[pair]:
                sp2 = listminus(sp[0], saved)
                if sp2 and sp2 not in akeys:
                    akeys.append(sp2)
        for sp in akeys:
            # kolla också för var och en i listan sp om de nollas av de andra i sp plus saved
            if len(sp) > 1:
                ok = True
                for p in sp:
                    ll = ciidic[p]
                    tset = listminus(sp,p) + saved
                    for ele,sig in ll:
                        if issubset(ele, tset):
                            ok = False
                            break
                if not ok:
                    continue
            saved = origsaved + sp
            remaining = listminus(origremaining, sp)
            removed = []
            changed = movetoremoved(remaining, saved, removed, ciidic)
            while changed:
                changed = False
                if movetosaved(remaining, saved, removed, ciidic):
                    changed = True
                if movetoremoved(remaining, saved, removed, ciidic):
                    changed = True
            if not remaining:
                s = set(saved)
                if s not in alternatives:
                    alternatives.append(s)
        if not alternatives:
            alternatives = [set(saved)]
            print("Failed to find clean condition alternatives")
            print("Saved: ", saved)
            print("Remaining: ", remaining)
            print("Dict: ", ciidic)
    else:
        alternatives = [set(saved)]
    resdic2lst = []
    for saved in alternatives:
        # preparera resdic2 från saved, dvs betinga var och en på övriga
        saved = list(saved)
        resdic2 = {}
        for pair in lst:
            if pair in saved:
                cond = allexcept(saved, pair)
            else:
                mx = ((),0.0)
                for sp in ciidic[pair]:
                    if issubset(sp[0], saved) and sp[1]>mx[1]:
                        mx = sp
                if mx[1] > 0.0:
                    cond = mx[0]
                else:
                    cond = allexcept(saved, pair)
            #(dic1, dic2, kdic) = get_cond_stats_one(df, pair, cond, selg, tags, eff)
            (dic1, dic2, kdic) = get_cond_stats_daydiff(df, pair, cond, selg, tags, eff)
            #tmp = analysediff(dic2, dic1, kdic)
            tmp = analysediff01(dic2, dic1, kdic)
            tmp['cond'] = cond
            resdic2[pair] = tmp
        for s in saved:
            if not s in totsaved:
                totsaved.append(s)
        resdic2lst.append(resdic2)
    return resdic2lst, totsaved

def analyse_one_effect(df, selg, pair, cond, eff):
    tags = []
    for tag in ["death", "inc", "nic"]:
        if eff in valsincolumn(df, 'event_time_' + tag):
            tags.append(tag)
    #(dic1, dic2, kdic) = get_cond_stats_one(df, pair, [], selg, tags, eff)
    (dic1, dic2, kdic) = get_cond_stats_daydiff(df, pair, [], selg, tags, eff)
    res0 = analysediff(dic2, dic1, kdic)
    #(dic1, dic2, kdic) = get_cond_stats_one(df, pair, cond, selg, tags, eff)
    (dic1, dic2, kdic) = get_cond_stats_daydiff(df, pair, cond, selg, tags, eff)
    res1 = analysediff(dic2, dic1, kdic)
    print("Correlates to " + eff)
    display_analysis_one_row(pair, res1, res0)
    print()
    #return { 'mean0': res0['mean'], 'sig0':res0['sig'], 'mean1': res1['mean'], 'sig1': res1['sig']}

def showeffects1(win, df, selg, coldic, tags, eff):
    name = "Correlates of " + str(list(coldic.keys())) + " to " + eff
    resdic = analyseeffects1(df, selg, coldic, tags, eff)
    show_analysis_rows_dict(win, name, resdic, False)

def show_all_effects(win, df, selg, sig, eff):
    coldic = {}
    selcols = [s[0] for s in selg if s[0] is not False]
    if 'sex' not in selcols:
        coldic['sex'] = ['flicka']  # special since using both are redundant
    for col in ['other_dia','diagnosis_class','surgery_diff','radio_diff','cytoclass_diff','stemcell_diff']:
        if col not in selcols:
            coldic[col] = valsincolumn(df, col)
    tags = []
    for tag in ["death", "inc", "nic"]:
        if eff in valsincolumn(df, 'event_time_' + tag):
            tags.append(tag)
    name = "Direct correlates to " + eff
    resdic1 = analyseeffects1(df, selg, coldic, tags, eff)
    resdic2lst,remaining = analyseeffects2new(df, selg, resdic1, sig, tags, eff)
    if len(resdic2lst) == 1:
        show_analysis_rows_dict(win, name, resdic2lst[0], resdic1)
    else:
        print("There are %d alternatives. Press return to switch." % len(resdic2lst))
        for i,resdic2 in enumerate(resdic2lst):
            show_analysis_rows_dict(win, name, resdic2lst[i], resdic1)
            if (i<len(resdic2lst)-1):
                input()

def display_all_effects(df, selg, sig, efflst = False, extratext="", tags = False):
    coldic = {}
    selcols = [s[0] for s in selg if s[0] is not False]
    if 'sex' not in selcols:
        coldic['sex'] = ['flicka']  # special since using both are redundant
    for col in ['other_dia','diagnosis_class','surgery_diff','radio_diff','cytoclass_diff','stemcell_diff']:
        if col not in selcols:
            coldic[col] = valsincolumn(df, col)
    effdic = {tag : valsincolumn(df, 'event_time_' + tag) for tag in (["death", "inc", "nic"] if tags is False else tags)}
    if efflst is False:
        efflst = []
        for tag in effdic:
            for eff in effdic[tag]:
                if eff not in efflst:
                    efflst.append(eff)
        efflst.sort()
    for eff in efflst:
        taglst = []
        for tag in effdic:
            if eff in effdic[tag]:
                taglst.append(tag)
        name = "Direct correlates to " + eff + extratext
        resdic1 = analyseeffects1(df, selg, coldic, taglst, eff)
        resdic2lst,remaining = analyseeffects2new(df, selg, resdic1, sig, taglst, eff)
        if len(resdic2lst) == 1:
            display_analysis_rows_dict(name, resdic2lst[0], resdic1)
            print()
        else:
            print("There are %d alternatives:" % len(resdic2lst))
            for i,resdic2 in enumerate(resdic2lst):
                display_analysis_rows_dict(name, resdic2lst[i], resdic1)
                print("--------------------")
            display_effect_all_comb(df, selg, remaining, eff)

def get_list_stats2(df, efftag, eff):
    # return min, max, first tsum, first nsum
    mn = 0
    mx = 0
    nsum_f = 0
    tsum_f = 0
    ncol = 'no_event_time_' + efftag
    ecol = 'event_time_' + efftag
    for i in range(len(df)):
        s = df.iloc[i]
        r = s[ncol]
        lst = s[ecol]
        ok = False
        for ele in lst:
            if (ele[0] == eff):
                if ele[1] < mn:
                    mn = ele[1]
                if ele[1] > mx:
                    mx = ele[1]
                tsum_f += ele[1]
                nsum_f += 1
                ok = True
        if not ok:
            tsum_f += r if not isnan(r) else 0
    return (mn, mx, tsum_f, nsum_f)

def show_profile_data(win, df, selg, pair, cond, eff):
    # först utan betingning
    (d1, d2) = select_data(df, [pair], False, selg)
    stats = {}
    for efftag in ['nic', 'inc', 'death']:
        stats1[(1,efftag)] = get_list_stats2(d1, efftag, eff)
        stats2[(2,efftag)] = get_list_stats2(d2, efftag, eff)
    mn = min([stats[k][0] for k in stats])
    mx = max([stats[k][1] for k in stats])
    win.fig.clear()
    show_profile_axes(win, 100, 200, 800, 100, mn, mx, 365)
    for (dind, dd, c1, c2, c3) in zip([False,True], [d1,d2], ['pink','lightgreen'], ['red','#00CD00'], ['#8B0000','#006400']):
        tsum = 0
        nsum = 0
        lst1 = []
        lst2 = []
        lst3 = []
        for i in range(len(dd)):
            s = dd.iloc[i]
            tacc = max([s['no_event_time_nic'], s['no_event_time_inc'], s['no_event_time_death']]) 
            nacc = 0
            for ele in s['event_time_nic']:
                if (ele[0] == eff):
                    lst1.append(ele[1])
                    if ele[1] < tacc:
                        tacc = ele[1]
                        nacc = 1
            for ele in s['event_time_inc']:
                if (ele[0] == eff):
                    lst2.append(ele[1])
                    if ele[1] < tacc:
                        tacc = ele[1]
                        nacc = 1
            for ele in s['event_time_death']:
                if (ele[0] == eff):
                    lst3.append(ele[1])
                    if ele[1] < tacc:
                        tacc = ele[1]
                        nacc = 1
            tsum += tacc
            nsum += nacc
        ft = tsum / nsum if nsum > 0 else None
        show_profile_bars(win, 100, 200, 800, 100, mn, mx, lst1, c1, c1, dind, 0)
        show_profile_bars(win, 100, 200, 800, 100, mn, mx, lst2, c2, c2, dind, 0)
        show_profile_bars(win, 100, 200, 800, 100, mn, mx, lst3, c3, c3, dind, 0)
        show_profile_bars(win, 100, 200, 800, 100, mn, mx, [ft], c2, c2, dind, 1)

    
def printsorteddicts(diclst):
    iset = set()
    for dic in diclst:
        iset = iset.union(dic.keys())
    ilst = sorted(list(iset))
    fstr = "%s" + ("\t%d" * len(diclst))
    for ele in ilst:
        tmp = [ele] + [dic[ele] if ele in dic else 0 for dic in diclst]
        print(fstr % tuple(tmp))

# Full genomlysning:
# relativt seneffekt eff och konkurerande grundorsaker pairlist för varje
# förekommande kombination av pair visa:
# hur många som fått eff av hur många som fått behandlingen och hur lång total obstid.
# visualisera individuella tider och medelfrekvensen

def display_effect_all_comb(df, selg, conds, eff, tags = False):
    effdic = {tag : valsincolumn(df, 'event_time_' + tag) for tag in (["death", "inc", "nic"] if tags is False else tags)}
    taglst = []
    for tag in effdic:
        if eff in effdic[tag]:
            taglst.append(tag)
    name = "Statistics of " + eff
    #(dic, datadic) = get_cond_stats_comb(df, conds, selg, taglst, eff)
    (dic, dummy1, dummy2) = get_cond_stats_daydiff(df, [], conds, selg, taglst, eff)
    print(name)
    tscale = 1.0/36525
    sum0 = 0
    sum1 = 0
    sum2 = 0
    for k in dic:
        sum0 += dic[k][0]
        sum1 += dic[k][1]
        sum2 += dic[k][2]
    pr1 = sum1/sum0
    pr2 = sum2/sum0
    for k in dic:
        if dic[k][2] > 0:
            resstr = ""
            for kk in k.split('_'):
                if kk:
                    if kk[0]=='~':
                        resstr += " "+kk
                    else:
                        resstr += "  "+kk
            resstr += "   "
            f1 = (dic[k][0] + 1)/((dic[k][1] + pr1)*tscale)
            f2 = (dic[k][0] + 1)*pr1/(dic[k][1] + pr1)
            resstr += "%c  %f  %.2f  ( %d / %d : %.2f )" % ('+' if f2 >= 1.5 else ' ', f1, f2, dic[k][0], dic[k][2], dic[k][1]*tscale)
            #if dic[k][1] > 0:
            #    resstr += "%f  ( %d / %d : %.2f )" % (dic[k][0]/(dic[k][1]*tscale),dic[k][0],dic[k][2],dic[k][1]*tscale)
            #else:
            #    resstr += "%f  ( %d / %d : %.2f )" % (0.0, dic[k][0],dic[k][2],dic[k][1]*tscale)
            print(resstr)
    n = 0
    for kk in k.split('_'):
        if kk:
            n += len(kk)+1
    resstr = " average " + " "*max(0, n-3)
    resstr += "%f  %.2f  ( %d / %d : %.2f )" % (1/(pr1*tscale), 1.0, 1, round(pr2), pr1*tscale)
    print(resstr)
    print()

def display_condensed_data(df, selg, eff):
    effcols = ['event_time_' + tag for tag in ["death", "inc", "nic"]]
    df = select_data_alt(df, [(col, eff) for col in effcols], selg)
    for i in range(len(df)):
        vec = df.iloc[i]
        resid = [c[0] for c in vec['event_time_inc'] if c[0][0]=='C']
        resstr = "%s\t%s\t%s" % (vec['LopNr'],vec['diagnosis'],("("+",".join(resid)+")") if resid else "")
        for col in ['cytoclass_diff','radio_diff','surgery_diff','stemcell_diff','other_dia']:
            for b in vec[col]:
                if type(b)==list:
                    resstr += "\t" + b[0]
                else:
                    resstr += "\t" + b
        print(resstr)

def collect_days(df, col):
    dic={}
    #cols = ['surgery_diff','cytoclass_diff','stemcell_diff','radio_diff']
    for lst in df[col]:
        for ele in lst:
            if not ele[0] in dic:
                dic[ele[0]] = (0,0)
            dic[ele[0]] = tupleadd(dic[ele[0]], (1,ele[1]))
    for key in dic:
        pair = dic[key]
        dic[key] = pair[1]/pair[0] if pair[0]>0 else 0
    return dic

def ddate_hist_dict(df, mn, mx):
    hdict = { k : 0 for k in range(mn, mx+1)}
    for i in range(len(df)):
        dd = df.iloc[i]['DiagnosDat']
        if dd in hdict:
            hdict[dd] += 1
        else:
            hdict[dd] = 1
    return hdict

def show_dd_axes(win, x, y, wdt, hgt, labels, count):
    win.add_line(x, y, x+wdt, y, 0.75, 'black')
    win.add_line(x, y+hgt+1, x+wdt, y+hgt+1, 0.75, 'black')
    for pr in [0.25, 0.5, 0.75]:
        ln = win.add_line(x, y+pr*hgt+1, x+wdt, y+pr*hgt+1, 0.75, 'lightgray')
        ln.set_zorder(1)
    onewdt = min(wdt/len(labels), hgt)
    for ind in range(len(labels)):
        win.add_text(x+onewdt/2+ind*onewdt, y-20, labels[ind], 12, 'center')
    win.add_text(x-8, y-4, "0", 9, 'right')
    win.add_text(x-8, y+hgt-4, "%d" % count, 9, 'right')

def show_dd_bar(win, x, y, wdt, hgt, num, ind, maxcount, val, col):
    onewdt = min(wdt/num, hgt)
    bwdt = onewdt-2
    boff = 1
    win.add_rect(x + onewdt*ind + boff, y, bwdt, hgt*val/maxcount, 0, None, col)

def show_dd_data(win, dic, x1, x2, mx, pos):
    win.fig.clear()
    show_histogram_axes(win, 60, 50+pos*250, 800, 200, ['1970','2015'], 100)
    for k in dic:
        show_dd_bar(win, 60, 50+pos*250, 800, 200, 46, int(k)-1970, 100, dic[k], 'orange')

    

def show_dd_data__(win, dic, sel1, sel3, selg, letters):
    dd = select_data3(df, sel1, sel3, selg)
    win.fig.clear()
    show_histogram_axes(win, 100, 100, 800, 200, letters, len(dd[0]))
    cols = ['#8B0000','red','pink',None]
    #alpha = 0.0
    #counts = [alpha/4, alpha/4, alpha/4, alpha/4]
    for ind in range(len(letters)):
        letter = letters[ind]
        counts = [0.0, 0.0, 0.0, 0.0]
        for i in range(len(dd[0])):
            s = dd[0].iloc[i]
            ok = False
            for ele in s['event_time_death']:
                if (letter is False or ele[0] == letter):
                    counts[0] += 1
                    ok = True
            if not ok:
                for ele in s['event_time_inc']:
                    if (letter is False or ele[0] == letter):
                        counts[1] += 1
                        ok = True
            if not ok:
                for ele in s['event_time_nic']:
                    if (letter is False or ele[0] == letter):
                        counts[2] += 1
                        ok = True
            if not ok:
                counts[3] += 1
        total = sum(counts)
        normcount = 0.0
        normtotal = 0.0
        for i in range(len(dd[2])):
            s = dd[2].iloc[i]
            ok = False
            for ele in s['event_time_death']:
                if (letter is False or ele[0] == letter):
                    normcount += 1
                    ok = True
            if not ok:
                for ele in s['event_time_inc']:
                    if (letter is False or ele[0] == letter):
                        normcount += 1
                        ok = True
            if not ok:
                for ele in s['event_time_nic']:
                    if (letter is False or ele[0] == letter):
                        normcount += 1
                        ok = True
            normtotal += 1
        show_histogram_bars(win, 100, 100, 800, 200, len(letters), ind, [x/total for x in counts], normcount/normtotal,cols)

def add_complex_attribute(df, col, name, oldlst):
    for i in range(len(df)):
        lst = df.iloc[i][col]
        x = [ ele[1] for ele in lst if ele[0] in oldlst]
        if x:
            lst.append([name, min(x)])
