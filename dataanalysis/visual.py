"""
Copyright (C) 2018-2021 RISE Research Institute of Sweden AB

File: visual.py

Author: anders.holst@ri.se

"""

from math import *
from hist import *
from window import *
from Color import *

def winresize(win, width, height):
    win.width = width
    win.height = height
    win.fig.set_size_inches((width/100.0, height/100.0))
    win.trans = win.fig.transFigure.inverted()

def calctickdist(span, pixdist, goaldist):
    onespan = span * goaldist / pixdist
    lsp = log(onespan)/log(10)+0.15
    fact = [1.0, 2.0, 5.0][int(floor((lsp - floor(lsp))*3))]
    return pow(10, floor(lsp))*fact
    
def show_hist(win, x, y, wdt, hgt, hist, hx1, hx2, col):
    ymax = max(hist[1]) * 1.05
    yscale = hgt/ymax if ymax > 0 else 0
    xscale = wdt/(hx2 - hx1)
    hlist = list(map(lambda xx,yy: (x + (xx-hx1)*xscale, y + yy*yscale), hist[0], hist[1]))
    hlist = [[x + (hist[0][0]-hx1)*xscale, y]] + hlist + [[x + (hist[0][-1]-hx1)*xscale, y]]
    win.add_polygon(hlist, 0, False, col)

def show_hist_axes(win, x, y, wdt, hgt, desc, hx1, hx2, tickdist):
    mx = floor(hx2/tickdist)
    mn = ceil(hx1/tickdist)
    xscale = wdt/(hx2 - hx1)
    win.add_line(x, y, x+wdt, y, 0.75, 'black')
    x0 = x - xscale*hx1
    for i in range(mn, mx+1):
        xx = x0 + i*tickdist*xscale
        win.add_line(xx, y-desc, xx, y+hgt if i==0 else y, 0.75, 'black')
        win.add_text(xx, y-desc-12, "%g" % (i*tickdist), 9, 'center')

def show_hist_axes_v2(win, x, y, wdt, hgt, desc, hx1, hx2, tickdist):
    mx = floor(hx2/tickdist)
    mn = ceil(hx1/tickdist)
    xscale = wdt/(hx2 - hx1)
    win.add_line(x, y, x+wdt, y, 0.75, 'black')
    x0 = x - xscale*hx1
    for i in range(mn, mx+1):
        xx = x0 + i*tickdist*xscale
        win.add_line(xx, y-desc, xx, y+hgt if i==0 else y, 0.75, 'black')
        win.add_text(xx, y-desc-18, "%.1f" % (i*tickdist), 12, 'center')

def show_hist_stats(win, x, y, wdt, hgt, hist, hx1, hx2, eps, col):
      xscale = wdt/(hx2 - hx1)
      mean, var = hist_mean_var(hist)
      xmin, xmax = hist_maxinterval(hist, 1.0-eps)
      win.add_line(x + (xmin-hx1)*xscale, y+int(hgt/3), x + (xmax-hx1)*xscale, y+int(hgt/3), 0.75, col)
      win.add_line(x + (xmin-hx1)*xscale, y+int(hgt/6), x + (xmin-hx1)*xscale + 0.01, y+int(hgt/2), 0.75, col)
      win.add_line(x + (xmax-hx1)*xscale, y+int(hgt/6), x + (xmax-hx1)*xscale + 0.01, y+int(hgt/2), 0.75, col)
      win.add_line(x + (mean-hx1)*xscale, y, x + (mean-hx1)*xscale+0.01, y+int(hgt*2/3), 0.75, col)

      
# 1) analysis row: Effect key, signif, diff hist, abs hist
# 2) analysis cond row: Effect key, signif, diff hist and cond, abs hist
# 3) analysis grid: signif and mean
# 4) analysis cond grid: signif and mean and cond
# Input: straight or cond stats

lowsign = 0.01
highsign = 0.0001

def set_sign_level(low, high):
    global lowsign
    global highsign
    lowsign = low
    highsign = high

def signcolor(sign):
    global lowsign
    global highsign
    if sign > lowsign:
        return 'white'
    elif sign < highsign:
        return 'black'
    else:
        llow = log(lowsign)
        lhigh = log(highsign)
        gr = (log(sign) - llow)/(lhigh - llow)
        return hsl_color(0.0, 0.0, 0.5-1.5*gr)

def show_analysis_rows(win, title, dic, dic0 = False):
    win.fig.clear()
    maxd = max([abs(dic[k]['mean']) + 2*dic[k]['std'] for k in dic])
    maxl = max([max(dic[k]['mean1'], dic[k]['mean2']) + 2*dic[k]['std'] for k in dic])
    if dic0 is not False:
        maxd = max(maxd, max([abs(dic0[k]['mean']) + 2*dic0[k]['std'] for k in dic0]))
        maxl = max(maxl, max([max(dic0[k]['mean1'], dic0[k]['mean2']) + 2*dic0[k]['std'] for k in dic0]))
    num = len(dic)
    rowhgt = int(win.height / (num + 1))
    histwdt1 = int((win.width - 240)*2*maxd/(2*maxd+maxl))
    histwdt2 = int((win.width - 240)*maxl/(2*maxd+maxl))
    tickdist = pow(10, floor(log(maxl/2)/log(10)))
    cold = hsl_color(0.35, 1.0, -0.1)
    cold0 = hsl_color(0.35, 1.0, 0.4)
    col1 = hsl_color(0.16, 1.0, 0.2)
    col2 = hsl_color(0.58, 1.0, 0.2)
    win.add_text(win.width/2, win.height - rowhgt/4 - 8, title, 12, 'center')
    show_hist_axes(win, 160, rowhgt/2, histwdt1, win.height - rowhgt, 6, -maxd, maxd, tickdist)
    show_hist_axes(win, 200+histwdt1, rowhgt/2, histwdt2, win.height - rowhgt, 6, 0.0, maxl, tickdist)
    ind = 0
    for k in sorted(dic):
        win.add_text(20, win.height - rowhgt*(1+ind) - 8, str(k), 12)
        if dic0 is not False:
            win.add_ellipse(58, win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic0[k]['sig']))
        win.add_ellipse(50, win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic[k]['sig']))
        win.add_text(80, win.height - rowhgt*(1+ind) - 6, "%.6f" % dic[k]['sig'], 10)
        if dic[k]['hist'] is not False:
            (mn, mx) = hist_maxinterval(dic[k]['hist'], 1.0-lowsign)
            pk = hist_peak(dic[k]['hist'])
            win.add_text(160 + histwdt1*(pk/maxd + 1.0)*0.5, win.height - rowhgt*(1+ind) - 6, "[%.2f,%.2f,%.2f]" % (mn, pk, mx), 10, 'center')
            if dic0 is not False:
                show_hist(win, 160, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic0[k]['hist'], -maxd, maxd, cold0)
            show_hist(win, 160, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic[k]['hist'], -maxd, maxd, cold)
        if dic[k]['hist1'] is not False:
            show_hist(win, 200+histwdt1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist1'], 0.0, maxl, col1)
            win.add_text(200+histwdt1+histwdt2*dic[k]['mean1']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n1'],dic[k]['t1']), 10, 'left' if dic[k]['mean1']>dic[k]['mean2'] else 'right')
        if dic[k]['hist2'] is not False:
            show_hist(win, 200+histwdt1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist2'], 0.0, maxl, col2)
            win.add_text(200+histwdt1+histwdt2*dic[k]['mean2']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n2'],dic[k]['t2']), 10, 'right' if dic[k]['mean1']>dic[k]['mean2'] else 'left')
        ind += 1

def show_analysis_rows_dict(win, title, dic, dic0 = False):
    win.fig.clear()
    num = len(dic)
    rowhgt = int(win.height / (num + 1))
    win.add_text(win.width/2, win.height - rowhgt/4 - 8, title, 12, 'center')
    if num == 0:
        return
    maxd = max([abs(dic[k]['mean']) + 2*dic[k]['std'] for k in dic])
    maxl = max([max(dic[k]['mean1'], dic[k]['mean2']) + 2*dic[k]['std'] for k in dic])
    if dic0 is not False:
        maxd = max(maxd, max([abs(dic0[k]['mean']) + 2*dic0[k]['std'] for k in dic0]))
        maxl = max(maxl, max([max(dic0[k]['mean1'], dic0[k]['mean2']) + 2*dic0[k]['std'] for k in dic0]))
    histwdt1 = int((win.width - 240)*2*maxd/(2*maxd+maxl))
    histwdt2 = int((win.width - 240)*maxl/(2*maxd+maxl))
    tickdist = pow(10, floor(log(maxl/2)/log(10)))
    cold = hsl_color(0.35, 1.0, -0.1)
    cold0 = hsl_color(0.35, 1.0, 0.4)
    col1 = hsl_color(0.16, 1.0, 0.2)
    col2 = hsl_color(0.58, 1.0, 0.2)
    show_hist_axes(win, 160, rowhgt/2, histwdt1, win.height - rowhgt, 6, -maxd, maxd, tickdist)
    show_hist_axes(win, 200+histwdt1, rowhgt/2, histwdt2, win.height - rowhgt, 6, 0.0, maxl, tickdist)
    ind = 0
    for k in sorted(dic):
        if dic0 is not False and k in dic0:
            win.add_ellipse(28, win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic0[k]['sig']))
        win.add_ellipse(20, win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic[k]['sig']))
        win.add_text(50, win.height - rowhgt*(1+ind) - 6, "%.4f" % dic[k]['sig'], 10)
        nm = str(k[1]) if type(k)==tuple else str(k)
        if 'cond' in dic[k]:
            k2 = dic[k]['cond']
            nm += " | " + (k2[1] if type(k2)==tuple else str(k2))
        win.add_text(120, win.height - rowhgt*(1+ind) - 8, nm, 12)
        if dic[k]['hist'] is not False:
            (mn, mx) = hist_maxinterval(dic[k]['hist'], 1.0-lowsign)
            pk = hist_peak(dic[k]['hist'])
            win.add_text(160 + histwdt1*(pk/maxd + 1.0)*0.5, win.height - rowhgt*(1+ind) - 6, "[%.2f,%.2f,%.2f]" % (mn, pk, mx), 10, 'center')
            if dic0 is not False and k in dic0:
                show_hist(win, 160, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic0[k]['hist'], -maxd, maxd, cold0)
            show_hist(win, 160, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic[k]['hist'], -maxd, maxd, cold)
        if dic[k]['hist1'] is not False:
            show_hist(win, 200+histwdt1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist1'], 0.0, maxl, col1)
            win.add_text(200+histwdt1+histwdt2*dic[k]['mean1']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n1'],dic[k]['t1']), 10, 'left' if dic[k]['mean1']>dic[k]['mean2'] else 'right')
        if dic[k]['hist2'] is not False:
            show_hist(win, 200+histwdt1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist2'], 0.0, maxl, col2)
            win.add_text(200+histwdt1+histwdt2*dic[k]['mean2']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n2'],dic[k]['t2']), 10, 'right' if dic[k]['mean1']>dic[k]['mean2'] else 'left')
        ind += 1

def show_analysis_rows_dict_v2(win, title, dic, dic0 = False, params = {}):
    paramdict = {'rowheight':50, 'namesize': 120, 'marg0': 0, 'marg1': 60, 'marg2' : 40, 'shownum':True, 'col1': 0.16, 'col2': 0.58, 'legendY': False, 'legendX1': False, 'legendX2': False }
    paramdict.update(params)
    win.fig.clear()
    num = len(dic)
    rowhgt = paramdict['rowheight'] # int(win.height / (num + 1))
    marg = paramdict['namesize'] + paramdict['marg0'] + paramdict['marg1'] + paramdict['marg2'] + 80
    if num == 0:
        return
    if 'maxd' in paramdict and 'maxl' in paramdict:
        maxd = paramdict['maxd']
        maxl = paramdict['maxl']
    else:
        maxd = max([abs(dic[k]['mean']) + 2*dic[k]['std'] for k in dic])
        maxl = max([max(dic[k]['mean1'], dic[k]['mean2']) + 2*dic[k]['std'] for k in dic])
        if dic0 is not False:
            maxd = max(maxd, max([abs(dic0[k]['mean']) + 2*dic0[k]['std'] for k in dic0]))
            maxl = max(maxl, max([max(dic0[k]['mean1'], dic0[k]['mean2']) + 2*dic0[k]['std'] for k in dic0]))
    if 'fixscale' in paramdict:
        newwidth = int(paramdict['fixscale']*(2*maxd+maxl) + marg)
        newheight = int(rowhgt*(num + 2.5))
        winresize(win, newwidth, newheight)
    histwdt1 = int((win.width - marg)*2*maxd/(2*maxd+maxl))
    histwdt2 = int((win.width - marg)*maxl/(2*maxd+maxl))
    tickdist = calctickdist(maxl, histwdt2, 80) #pow(10, floor(log(maxl/2)/log(10)))
    lpos = 20 + paramdict['marg0']
    bpos = lpos + paramdict['namesize']
    hpos1 = bpos + 20 + paramdict['marg1']
    hpos2 = hpos1 + paramdict['marg2'] + histwdt2
    cold = hsl_color(0.35, 1.0, -0.1)
    cold0 = hsl_color(0.35, 1.0, 0.4)
    col1 = hsl_color(paramdict['col1'], 1.0, 0.2)
    col2 = hsl_color(paramdict['col2'], 1.0, 0.2)
    win.add_text(win.width/2, win.height - rowhgt/4 - 8, title, 12, 'center')
    show_hist_axes_v2(win, hpos2, 2*rowhgt, histwdt1, win.height - rowhgt*2.5, 6, -maxd, maxd, tickdist)
    show_hist_axes_v2(win, hpos1, 2*rowhgt, histwdt2, win.height - rowhgt*2.5, 6, 0.0, maxl, tickdist)
    if paramdict['legendY']:
        txt_y = win.add_text(36, win.height/2, paramdict['legendY'], 12, 'center')
        txt_y.set_rotation_mode('anchor')
        txt_y.set_rotation('vertical')
    if paramdict['legendX1']:
        txt_x1 = win.add_text(hpos1, 12, paramdict['legendX1'], 12, 'left')
    if paramdict['legendX2']:
        txt_x2 = win.add_text(hpos2+histwdt1/2, 12, paramdict['legendX2'], 12, 'center')
    ind = 0
    for k in sorted(dic):
        if dic0 is not False and k in dic0:
            win.add_ellipse(28+paramdict['namesize'], win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic0[k]['sig']))
        win.add_ellipse(bpos, win.height - rowhgt*(1+ind) - 8, 16, 16, 1, 'gray', signcolor(dic[k]['sig']))
        #win.add_text(50, win.height - rowhgt*(1+ind) - 6, "%.4f" % dic[k]['sig'], 10)
        nm = str(k[1]) if type(k)==tuple else str(k)
        #if 'cond' in dic[k]:
        #    k2 = dic[k]['cond']
        #    nm += " | " + (k2[1] if type(k2)==tuple else str(k2))
        win.add_text(lpos, win.height - rowhgt*(1+ind) - 8, nm, 12)
        if dic[k]['hist'] is not False:
            (mn, mx) = hist_maxinterval(dic[k]['hist'], 1.0-lowsign)
            pk = hist_peak(dic[k]['hist'])
            if paramdict['shownum']:
                win.add_text(hpos2 + histwdt1*(pk/maxd + 1.0)*0.5, win.height - rowhgt*(1+ind) - 6, "[%.2f,%.2f,%.2f]" % (mn, pk, mx), 10, 'center')
            if dic0 is not False and k in dic0:
                show_hist(win, hpos2, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic0[k]['hist'], -maxd, maxd, cold0)
            show_hist(win, hpos2, win.height - rowhgt*(1.45+ind), histwdt1, rowhgt*0.9, dic[k]['hist'], -maxd, maxd, cold)
        if dic[k]['hist1'] is not False:
            show_hist(win, hpos1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist1'], 0.0, maxl, col1)
            if paramdict['shownum']:
                win.add_text(hpos1 + histwdt2*dic[k]['mean1']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n1'],dic[k]['t1']), 10, 'left' if dic[k]['mean1']>dic[k]['mean2'] else 'right')
        elif dic0 is not False and k in dic0 and dic0[k]['hist1'] is not False:
            show_hist(win, hpos1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic0[k]['hist1'], 0.0, maxl, col1)
            if paramdict['shownum']:
                win.add_text(hpos1 + histwdt2*dic0[k]['mean1']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic0[k]['n1'],dic0[k]['t1']), 10, 'left' if dic0[k]['mean1']>dic0[k]['mean2'] else 'right')
        if dic[k]['hist2'] is not False:
            show_hist(win, hpos1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic[k]['hist2'], 0.0, maxl, col2)
            if paramdict['shownum']:
                win.add_text(hpos1 + histwdt2*dic[k]['mean2']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic[k]['n2'],dic[k]['t2']), 10, 'right' if dic[k]['mean1']>dic[k]['mean2'] else 'left')
        elif dic0 is not False and k in dic0 and dic0[k]['hist2'] is not False:
            show_hist(win, hpos1, win.height - rowhgt*(1.45+ind), histwdt2, rowhgt*0.9, dic0[k]['hist2'], 0.0, maxl, col2)
            if paramdict['shownum']:
                win.add_text(hpos1 + histwdt2*dic0[k]['mean2']/maxl, win.height - rowhgt*(1+ind) - 6, "(%d/%.0f)" % (dic0[k]['n2'],dic0[k]['t2']), 10, 'right' if dic0[k]['mean1']>dic0[k]['mean2'] else 'left')
        ind += 1

def display_analysis_one_row(pair, res, res0 = False):
    resstr = ""
    nm = str(pair[1]) if type(pair)==tuple else str(pair)
    if 'cond' in res:
        k2 = res['cond']
        if type(k2) == list:
            nm += " | " + ",".join([kk[1] if type(kk)==tuple else str(kk) for kk in k2])
        elif type(k2) == tuple:
            nm += " | " + k2[1]
        else:
            nm += " | " + str(k2)
    sig = res['sig']
    mean = res['mean']
    fact = res['fact']
    fact00 = res['mean2']/res['mean1']
    #if res['hist'] is not False:
    #    pk = hist_peak(res['hist'])
    if sig < 0.0001:
        resstr += "***"
    elif sig < 0.001:
        resstr += "**"
    elif sig < 0.01:
        resstr += "*"
    else:
        resstr += ""
    resstr += "\t%.4f\t%.4f\t%.2f" % (sig, mean, fact)
    if res0 is not False:
        sig0 = res0['sig']
        mean0 = res0['mean']
        fact0 = res0['fact']
        resstr += "\t%.4f\t%.4f\t%.2f" % (sig0, mean0, fact0)
    else:
        resstr += "\t\t\t"
    resstr += "\t%s" % (nm)
    #sp = max(2, 50-len(resstr))
    resstr += "\t(%d/%d:%.0f)\t(%d/%d:%.0f)" % (res['n1'],res['s1'],res['t1'],res['n2'],res['s2'],res['t2'])
    print(resstr)

def display_analysis_rows_dict(title, dic, dic0 = False):
    num = len(dic)
    print(title)
    if num == 0:
        return
    for k in sorted(dic):
        display_analysis_one_row(k, dic[k], dic0[k] if dic0 is not False and k in dic0 else False)

def show_analysis_grid(win, title, gdic, gdic0 = False):
    win.fig.clear()
    numc = len(gdic)
    cols = sorted(gdic)
    numr = len(gdic[cols[0]])
    rows = sorted(gdic[cols[0]])
    rowhgt = int(win.height / (numr + 2))
    colwdt = int((win.width - 50)/ (numc if numc > 0 else 1))  
    win.add_text(win.width/2, win.height - rowhgt/4 - 8, title, 12, 'center')
    indc = 0
    for g in cols:
        win.add_text(50+indc*colwdt, win.height - rowhgt - 8, str(g), 12)
        indc +=1
    indr = 0
    for k in rows:
        win.add_text(20, win.height - rowhgt*(2+indr) - 8, str(k), 12)
        indc = 0
        for g in cols:
            if gdic0 is not False:
                win.add_ellipse(58+indc*colwdt, win.height - rowhgt*(2+indr) - 8, 16, 16, 1, 'gray', signcolor(gdic0[g][k]['sig']))
            win.add_ellipse(50+indc*colwdt, win.height - rowhgt*(2+indr) - 8, 16, 16, 1, 'gray', signcolor(gdic[g][k]['sig']))
            indc += 1
        indr += 1

def show_profile_axes(win, x, y, wdt, hgt, hx1, hx2, tickdist):
    mx = floor(hx2/tickdist)
    mn = ceil(hx1/tickdist)
    xscale = wdt/(hx2 - hx1)
    win.add_line(x, y, x+wdt, y, 0.75, 'black')
    x0 = x - xscale*hx1
    for i in range(mn, mx+1):
        xx = x0 + i*tickdist*xscale
        if i==0:
            win.add_line(xx, y-hgt, xx, y+hgt, 0.75, 'black')
        else:
            win.add_line(xx, y-4, xx, y+4, 0.75, 'black')
#        win.add_text(xx, y-desc-12, "%g" % (i*tickdist), 9, 'center')

def show_profile_bars(win, x, y, wdt, hgt, hx1, hx2, lst, col1, col2, below, full):
    xscale = wdt/(hx2 - hx1)
    first = True
    lst = [x for x in lst if x is not None]
    for val in lst:
        xx = x + xscale*(val - hx1)
        if first and val > 0:
            col = col1
            first = False
        else:
            col = col2
        ihgt = hgt*2/3 if not full else hgt
        if below:
            win.add_line(xx, y-4, xx, y-ihgt+4, 0.75, col)
        else:
            win.add_line(xx, y+4, xx, y+ihgt-4, 0.75, col)

def show_km_axes(win, x, y, wdt, hgt, hx1, hx2, xtickdist, ymx, ytickdist):
    mx = floor(hx2/xtickdist)
    mn = ceil(hx1/xtickdist)
    xscale = wdt/(hx2 - hx1)
    yscale = hgt/ymx
    win.add_line(x-4, y, x+wdt, y, 1, 'black')
    win.add_text(x-6, y-6, str(0), 12, 'right')
    for i in range(1,floor(ymx/ytickdist)+1):
        yy = y+yscale*ytickdist*i
        win.add_line(x-4, yy, x+wdt, yy, 0.75, 'lightgray')
        win.add_text(x-6, yy-6, str(int(ytickdist*i*100)), 12, 'right')
    x0 = x - xscale*hx1
    for i in range(mn, mx+1):
        xx = x0 + i*xtickdist*xscale
        if i==0:
            win.add_line(xx, y-4, xx, y+hgt+20, 1, 'black')
        else:
            win.add_line(xx, y-4, xx, y+4, 1, 'black')
        win.add_text(xx, y-20, str(i*xtickdist), 12, 'center')

def show_km_data(win, x, y, wdt, hgt, hx1, hx2, ymx, off, lst, col):
    xscale = wdt/(hx2 - hx1)
    yscale = hgt/ymx
    x0 = x - xscale*hx1 + xscale*off
    points = [(x0 + xscale*i, y + yscale*yy) for i,yy in enumerate(lst)]
    obj = win.add_polygon(points, 1, col, None)
    obj.set_closed(False)


    
