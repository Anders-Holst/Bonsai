"""
Copyright (C) 2018-2021 RISE Research Institute of Sweden AB

File: hist.py

Author: anders.holst@ri.se

"""

from math import *

def seconddegree(a, b, c, fmin, fmax):
    # solves a x^2 + b x + c = 0, and returns the closest real value in the interval [fmin, fmax]
    if a==0.0:
        f = -c/b
        return max(min(f, fmax), fmin)
    else:
        f0 = -b/(2*a)
        fd2 = f0*f0 - c/a
        if fd2<0:
            f = f0
        elif f0 < (fmin+fmax)/2:
            f = f0 + sqrt(fd2)
        else:
            f = f0 - sqrt(fd2)
        return max(min(f, fmax), fmin)

def hist_peak(hist):
    mx = 0.0
    mxind = False
    for i in range(len(hist[0])):
        if mx < hist[1][i]:
            mx = hist[1][i]
            mxind = hist[0][i]
    return mxind

def hist_significance(hist, val, peak):
    sum = 0.0
    if val > peak:
        for i in range(len(hist[0])):
            if hist[0][i] >= val:
                sum += hist[1][i]
    else:
        for i in range(len(hist[0])):
            if hist[0][i] <= val:
                sum += hist[1][i]
    return sum

def hist_mean_var(hist):
    sum = 0
    sqsum = 0
    for i in range(len(hist[0])):
        sum += hist[0][i]*hist[1][i]
        sqsum += hist[0][i]*hist[0][i]*hist[1][i]
    return (sum, sqsum-sum*sum)

def hist_maxinterval(hist, mass):
    i = 0
    j = len(hist[0])-1
    while i < j:
        a1 = hist[1][i-1] if i>0 else 0.0
        a2 = hist[1][i]
        b1 = hist[1][j+1] if j<len(hist[0])-1 else 0.0
        b2 = hist[1][j]
        if a2 < b2:
            ma = (a2+a1)/2
            mb = (a2*a2 - b1*b1)/(2*(b2-b1)) if b2>b1 else 0.0
            if mass + ma + mb > 1.0:
                break
            mass += ma
            i += 1
        else:
            mb = (b2+b1)/2
            ma = (b2*b2 - a1*a1)/(2*(a2-a1)) if a2>a1 else 0.0
            if mass + ma + mb > 1.0:
                break
            mass += mb
            j -= 1
    if i == 0:
        if j == len(hist[0])-1:
            return (hist[0][i], hist[0][j])
        else:
            b1 = hist[1][j+1]
            b2 = hist[1][j]
            delta = (b1+b2)/2 + mass - 1.0
            xb = (seconddegree(1.0, 0.0, -2*delta*(b1-b2) -b2*b2, b1, b2) - b2)/(b1 - b2)
            return (hist[0][i], hist[0][j]*(1.0-xb) + hist[0][j+1]*xb)
    else:            
        if j == len(hist[0])-1:
            a1 = hist[1][i-1]
            a2 = hist[1][i]
            delta = (a1+a2)/2 + mass - 1.0
            xa = (seconddegree(1.0, 0.0, -2*delta*(a1-a2) -a2*a2, a1, a2) - a2)/(a1 - a2)
            return (hist[0][i]*(1.0-xa) + hist[0][i-1]*xa, hist[0][j])
        else:
            a1 = hist[1][i-1]
            a2 = hist[1][i]
            b1 = hist[1][j+1]
            b2 = hist[1][j]
            delta = (a1+a2+b1+b2)/2 + mass - 1.0
            f = seconddegree(b1-b2+a1-a2, 0.0, -a2*a2*(b1-b2) - b2*b2*(a1-a2) - 2*(a1-a2)*(b1-b2)*delta, max(a1, b1), min(a2, b2))
            xa = (f - a2)/(a1 - a2)
            xb = (f - b2)/(b1 - b2)
            return (hist[0][i]*(1.0-xa) + hist[0][i-1]*xa, hist[0][j]*(1.0-xb) + hist[0][j+1]*xb)

