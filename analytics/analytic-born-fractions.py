#! /usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0,'../../../PanScales/scripts/')
from hfile import *
from scipy import integrate

pdf_uncert_hess = True
npdf_hess =  40
npdf_mc   = 100

use_tabulated_lo = True

class TabulatedLO(object):
    def __init__(self):
        self.dsigma_dpt_table = get_array("lo.res", "muf_1.0-ipdf_0")
        self.lptmin = np.log(self.dsigma_dpt_table[ 0,0])
        self.lptmax = np.log(self.dsigma_dpt_table[-1,0])
        self.npt    = self.dsigma_dpt_table.shape[0]

    def set_muf(self, muf):
        self.dsigma_dpt_table = get_array("lo.res", f"muf_{muf:.1f}-ipdf_0")

    def set_ipdf(self, ipdf):
        self.dsigma_dpt_table = get_array("lo.res", f"muf_1.0-ipdf_{ipdf}")

    def dsigma_dpt(self, pt, ymin=0.0, ymax=1.7, include_q=True, include_g=True):
        if pt<700:   return 0.0
        if pt>=6800: return 0.0
        lpt = np.log(pt)
        #print (lpt, self.lptmin, self.lptmax)
        xpt = (lpt-self.lptmin)/(self.lptmax-self.lptmin)*(self.npt-1)
        ipt = int(xpt)
        fpt = xpt-ipt
        #print (xpt, ipt, fpt)
        lo = 0.0
        hi = 0.0
        if include_q:
            lo += self.dsigma_dpt_table[ipt  ,1]
            hi += self.dsigma_dpt_table[ipt+1,1]
        if include_g:
            lo += self.dsigma_dpt_table[ipt  ,2]
            hi += self.dsigma_dpt_table[ipt+1,2]
        if lo<=0 or hi<=0: return (1-fpt)*lo+fpt*hi
        return np.exp((1-fpt)*np.log(lo)+fpt*np.log(hi))

if use_tabulated_lo:
    LO = TabulatedLO()
else:
    import lo
    LO = lo.LO(13600.0)

#--------------------------------------------------------------------------------
# analytic results
#
# see ./analytics for some extra details
#--------------------------------------------------------------------------------

CF=4.0/3.0
CA=3.0
nf=5.0
TR=0.5
Tf=nf*TR

zmin_def = 0.1
zmax_def = 0.2

# zmax = zmin+x, lim_(x\to 0)
# 
# Iqq = CF*(2*np.log((1-zmin)/(1-zmax)) - (zmax-zmin)*(1+(zmax+zmin)/2))
#     = x * CF*(1+z^2)/(1-z)
# Igq = CF*(2*np.log(zmax/zmin) - 2*(zmax-zmin)*(1-(zmax+zmin)/4))
#     = x * CF*(1+(1-z)^2)/z
# 
# Iqg = 2*Tf/3.0 * ((zmax**3-zmin**3) + ((1-zmin)**3-(1-zmax)**3))
#     = x * 2*Tf * (z^2 + (1-z)^2)
# Igg = 2*CA*(np.log(zmax*(1-zmin)/(zmin*(1-zmax))) - 2*(zmax-zmin) + (zmax**2-zmin**2)/2 - (zmax**3-zmin**3)/3)
#     = x* 2*CA*(z/(1-z) + (1-z)/z + z*(1-z))


# DFLAP splitting functions
Pqq = lambda z: CF*(1+z*z)/(1-z)
Pgq = lambda z: CF*(1+(1-z)*(1-z))/z
Pqg = lambda z: 2*Tf * (np.square(z) + np.square(1-z))
Pgg = lambda z: 2*CA * ((1-z)/z + z/(1-z) + z*(1-z))

# integrals with fixed parent pt and range of z fraction
#
# Some of the expressions below assume zmin<zmax!!
Iqq = lambda zmin, zmax: CF*(2*np.log((1-zmin)/(1-zmax)) - (zmax-zmin)*(1+(zmax+zmin)/2))
Igq = lambda zmin, zmax: CF*(2*np.log(zmax/zmin) - 2*(zmax-zmin)*(1-(zmax+zmin)/4))
Iqg = lambda zmin, zmax: 2*Tf/3.0 * ((zmax**3-zmin**3) + ((1-zmin)**3-(1-zmax)**3))
Igg = lambda zmin, zmax: 2*CA*(np.log(zmax*(1-zmin)/(zmin*(1-zmax))) - 2*(zmax-zmin) + (zmax**2-zmin**2)/2 - (zmax**3-zmin**3)/3)

ffg_q = lambda zmin, zmax: Igq(zmin, zmax)/(Iqq(zmin, zmax)+Igq(zmin, zmax)) if zmin<zmax else Pgq(zmin)/(Pqq(zmin)+Pgq(zmin))
ffg_g = lambda zmin, zmax: Igg(zmin, zmax)/(Iqg(zmin, zmax)+Igg(zmin, zmax)) if zmin<zmax else Pgg(zmin)/(Pqg(zmin)+Pgg(zmin))

fg_q_groom = np.vectorize(ffg_q)
fg_g_groom = np.vectorize(ffg_g)

# integrals with fixed pt cuts
#
# Leading jet has at least "ptcut"
# Subleading jet is between ptmin and ptmax
#
# first a helper
def Jab_int(pt, Iab, include_q, include_g, ptcut, ptmin, ptmax):
    zmin = ptmin/pt
    zmax = min(ptmax/pt, 1-ptcut/pt)
    if zmin>=zmax: return 0.0
    res= Iab(zmin,zmax)*LO.dsigma_dpt(pt, ymin=0.0, ymax=1.7,
                                      include_q=include_q, include_g=include_g)
    #print (f"{pt=}, {res=}")
    return res

def Jab(Iab, include_q, include_g, ptcut, ptmin, ptmax):
    epsrel = 1e-6 if use_tabulated_lo else 1e-4
    epsabs = 1e-6 if use_tabulated_lo else 1e-4
    res,err= integrate.quad(Jab_int, ptcut+ptmin, 6800.0,
                            epsrel=epsrel, epsabs=epsabs,
                            args=(Iab, include_q, include_g, ptcut, ptmin, ptmax))
    return res

def ffg_q_dijet(ptcut, ptmin, ptmax):
    Jgq = Jab(Igq, True, False, ptcut, ptmin, ptmax)
    Jqq = Jab(Iqq, True, False, ptcut, ptmin, ptmax)
    return Jgq/(Jqq+Jgq)

def ffg_g_dijet(ptcut, ptmin, ptmax):
    Jgg = Jab(Igg, False, True, ptcut, ptmin, ptmax)
    Jqg = Jab(Iqg, False, True, ptcut, ptmin, ptmax)
    return Jgg/(Jqg+Jgg)

fg_q_dijet = np.vectorize(ffg_q_dijet)
fg_g_dijet = np.vectorize(ffg_g_dijet)

def f_with_err(f, ptcut, ptmin, ptmax):
    LO.set_muf(1.0)
    central = f(ptcut, ptmin, ptmax)
    res = np.zeros((central.shape[0],7))
    for i in range(5): res[:,i] = central

    # factorisation scale variation
    for muf in [0.5,2.0]:
        LO.set_muf(muf)
        r = f(ptcut, ptmin, ptmax)
        res[:,3] = np.minimum(res[:,3],r)
        res[:,4] = np.maximum(res[:,4],r)

    res[:,3] = res[:,3] - central
    res[:,4] = res[:,4] - central

    # factorisation scale variation
    for ipdf in range(1,npdf_hess):
        LO.set_ipdf(ipdf)
        r = f(ptcut, ptmin, ptmax)
        res[:,6] += np.square(r-central)
    res[:,6] = np.sqrt(res[:,6])
    res[:,5] = -res[:,6]

    # full variation
    res[:,1]=-np.sqrt(np.square(res[:,3])+np.square(res[:,5]))
    res[:,2]=+np.sqrt(np.square(res[:,4])+np.square(res[:,6]))
        
    return res
#print (fg_q_dijet(750.0, 150.0, np.asarray([200.0,220.0])))
#sys.exit(0)




#--------------------------------------------------------------------------------
# tabulation
#--------------------------------------------------------------------------------
fout=open("analytic-born-fractions.res", "w")

# dijet selection
#
# We try two versions: one which uses a zfraction computed assuming
# the leading jet has exactly the min allowed value, and one which
# integrates over the LO spectrum.
all_vars = ["ptmin", "ptmax", "drmin"]

for v in all_vars:
    print (f"dijet, {v=}")
    ptlead=700.0
    zmin = 150.0/(ptlead+150.0)
    zmax = 200.0/(ptlead+200.0)
    if v=="ptmin":
        an_x = np.linspace(100.0, 199.9, 51)
        an_q = f_with_err(fg_q_dijet,ptlead, an_x, 200.0)
        an_g = f_with_err(fg_g_dijet,ptlead, an_x, 200.0)
        an0_q = fg_q_groom(an_x/(an_x+700.0),zmax)
        an0_g = fg_g_groom(an_x/(an_x+700.0),zmax)
    elif v=="ptmax":
        an_x = np.linspace(150.1, 300.0, 51)
        an_q = f_with_err(fg_q_dijet,ptlead, 150.0, an_x)
        an_g = f_with_err(fg_g_dijet,ptlead, 150.0, an_x)
        an0_q = fg_q_groom(zmin,an_x/(an_x+700.0))
        an0_g = fg_g_groom(zmin,an_x/(an_x+700.0))
    else:
        an_x = np.linspace(0.4, 1.2, 51)
        fq = f_with_err(fg_q_dijet,np.asarray([ptlead]), 150.0, 200.0)
        fg = f_with_err(fg_g_dijet,np.asarray([ptlead]), 150.0, 200.0)
        an_q = np.zeros((51,7))
        an_g = np.zeros((51,7))
        for i in range(51):
            an_q[i,:] = fq
            an_g[i,:] = fg
        an0_q = 0.0*an_x + fg_q_groom(zmin,zmax)
        an0_g = 0.0*an_x + fg_g_groom(zmin,zmax)

    print (f"# sel_dijet-dep_{v}", file=fout)
    print ( "# errors are given as q-,g-,q+,g+", file=fout)
    print (f"# columns: param, gfrac|q, gfrac|g,errs_tot,muf_vars,pdf_vars,gfrac|q(lead=700), gfrac|g(lead=700)", file=fout)
    for i in range(an_x.shape[0]):
        print (an_x[i], end=' ', file=fout)
        for j in range(7): print (an_q[i,j], an_g[i,j], end=' ', file=fout)
        print (an0_q[i], an0_g[i], file=fout)
    print ("", file=fout)
    print ("", file=fout)
    fout.flush()

# groomed selection
all_vars = ["zmin", "zmax", "drmin"]

for v in all_vars:
    print (f"groom, {v=}")

    if v=="zmin":
        an_x = np.linspace(0.05, 0.2, 51)
        an_q = fg_q_groom(an_x,0.2)
        an_g = fg_g_groom(an_x,0.2)
    elif v=="zmax":
        an_x = np.linspace(0.1, 0.5, 51)
        an_q = fg_q_groom(0.1,an_x)
        an_g = fg_g_groom(0.1,an_x)
    else:
        an_x = np.linspace(0.0, 1.2, 51)
        an_q = 0.0*an_x+fg_q_groom(0.1,0.2)
        an_g = 0.0*an_x+fg_g_groom(0.1,0.2)
        
    print (f"# sel_groom-dep_{v}", file=fout)
    print (f"# columns: param, gfrac|q, gfrac|g", file=fout)
        
    for i in range(an_x.shape[0]):
        print (an_x[i], an_q[i], an_g[i], file=fout)
    print ("", file=fout)
    print ("", file=fout)
    fout.flush()
