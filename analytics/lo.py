#!/usr/bin/env python3

import lhapdf
from me import *
from scipy import integrate

from math import cosh, exp, sqrt, pi, log

#-------------------------------------
class LO(object):
    # ctor
    def __init__(self, sqrts, mur=1.0, muf=1.0, gg2qq=True, gg2gg=True, qq2gg=True, qq2qq=True, qg2qg=True, use_ht=True, ipdf=0):
        self.sqrts = sqrts
        
        self.muf   = muf
        self.mur   = mur
        self.gg2qq = gg2qq
        self.gg2gg = gg2gg
        self.qq2gg = qq2gg
        self.qq2qq = qq2qq
        self.qg2qg = qg2qg
        
        self.use_ht = use_ht

        self.pdfs = lhapdf.mkPDFs("PDF4LHC21_40")
        self.ipdf = ipdf

        self.beam1_proton = True
        self.beam2_proton = True
        self.GeV2nb = 0.389e6

        self.epsabs = 1.0e-6
        self.epsrel = 1.0e-6
        self.xmax = 0.999
        
    # dsigma/(dpt dy) in nb/GeV: integrand (over Y=(y2-y1)/2)
    # note that a factor 16 \pi^2 \alphas^2 is taken out of the amplitudes and put in here
    def integrand_lo_dsigma_dpt_dy(self, Y, pt, y ,include_q, include_g):
        B1 = lambda i: i if self.beam1_proton else -i
        B2 = lambda i: i if self.beam2_proton else -i
        
        yhat = y-Y
        ch = cosh(yhat)

        x1 = 2.0*pt/self.sqrts*ch*exp(+Y)
        x2 = 2.0*pt/self.sqrts*ch*exp(-Y)

        if x1>=self.xmax: return 0.0
        if x2>=self.xmax: return 0.0

        t = -pt*pt*(1+exp(-2*yhat))
        u = -pt*pt*(1+exp(+2*yhat))
        s = -t-u

        scale = self.muf * pt;
        if self.use_ht: scale*=2
        xpdf1 = self.pdfs[self.ipdf].xfxQ(x1, scale)
        xpdf2 = self.pdfs[self.ipdf].xfxQ(x2, scale)

        #print (xpdf1)
        #print (xpdf2)

        nf=5
        
        sum=0.0

        # qq' (qqbar', qbar q' and qbarqbar') contrib
        pdf=0.0
        for i in range(1,nf+1):
            q1 = xpdf1[i]+xpdf1[-i]
            for j in range(1,nf+1):
                if i==j: continue
                pdf += q1*(xpdf2[j]+xpdf2[-j])
        if include_q:
            sum += pdf * (Mqqp2qqp(s,t,u) + Mqqp2qqp(s,u,t))

        # qq (and qbarqbar) contrib
        pdf=0.0
        for i in range(1,nf+1):
            pdf += xpdf1[B1(i)]*xpdf2[B2(i)] + xpdf1[B1(-i)]*xpdf2[B2(-i)]
        if include_q:
            sum += pdf * 2.0 * Mqq2qq(s,t,u)
 
        # qqbar (and qbarq) contrib
        pdf=0.0
        for i in range(1,nf+1):
            pdf += xpdf1[B1(i)]*xpdf2[B2(-i)] + xpdf1[B1(-i)]*xpdf2[B2(i)]
        if include_q:
            sum += pdf * (Mqqbar2qpqbarp(s,t,u) + Mqqbar2qpqbarp(s,u,t) + Mqqbar2qqbar(s,t,u) + Mqqbar2qqbar(s,u,t))
        if include_g:
            sum += pdf * 2 * Mqqbar2gg(s,t,u)

        # qg (gq, qbqrg and gqbqr) contrib
        pdf=0.0
        for i in range(1,nf+1):
            pdf += (xpdf1[i]+xpdf1[-i])*xpdf2[21]
        if include_q:
            sum += pdf * Mqg2qg(s,t,u)
        if include_g:
            sum += pdf * Mqg2qg(s,u,t)

        pdf=0.0
        for i in range(1,nf+1):
            pdf += xpdf1[21]*(xpdf2[-i]+xpdf2[i])
        if include_q:
            sum += pdf * Mqg2qg(s,u,t)
        if include_g:
            sum += pdf * Mqg2qg(s,t,u)
    
        # gg contrib
        pdf = xpdf1[21] * xpdf2[21]
        if include_g:
            sum += pdf*2.0*Mgg2gg(s,t,u)
        if include_q:
            sum += pdf*(Mgg2qqbar(s,t,u) + Mgg2qqbar(s,u,t))

        return sum/ch**4

    # spectrum as a function of py and y
    def dsigma_dpt_dy(self, pt, y, include_q=True, include_g=True):
        # impose boundaries on y
        tmp = self.sqrts/2.0/pt
        maxrap = log(tmp + sqrt(tmp*tmp-1))
        if abs(y) > maxrap: return 0.0

        Ymin = y-0.5*log(self.sqrts*exp(+y)/pt-1)
        Ymax = y+0.5*log(self.sqrts*exp(-y)/pt-1)

        res, err = integrate.quad(self.integrand_lo_dsigma_dpt_dy, Ymin, Ymax, args=(pt,y,include_q, include_g), epsabs=self.epsabs, epsrel=self.epsrel)

        scale = self.mur * pt
        if self.use_ht: scale*=2
        
        alphas = self.pdfs[self.ipdf].alphasQ(scale)
        return 0.25*pi*alphas**2*res/pt**3 * self.GeV2nb


    #-------------------------------------
    #dsigma/dpt
    def dsigma_dpt(self, pt, ymin=-1000, ymax=1000, include_q=True, include_g=True):
        if pt>=0.999*0.5*self.sqrts: return 0.0

        tmp = self.sqrts/2.0/pt
        maxrap = log(tmp + sqrt(tmp*tmp-1))
        ylo = -maxrap if ((ymax<ymin) or (ymin<-maxrap)) else ymin
        yup =  maxrap if ((ymax<ymin) or (ymax>+maxrap)) else ymax

        f = lambda y: self.dsigma_dpt_dy(pt, y, include_q, include_g)
        res, err = integrate.quad(f, ylo, yup)

        return res

# when run as a script, build a tabulation of the spectrum used in the
# dijet selection
if __name__ == "__main__":
    import sys
    import numpy as np
    
    sqrts = 13600.0
    mur   =     1.0

    print ("# columns are pt dsigma_dpt(q)  dsigma_dpt(g)")
    for muf in [1.0,0.5,2.0]:
        for ipdf in range (41):
            print (f"# muf_{muf}-ipdf_{ipdf}", file=sys.stderr)
            print (f"# muf_{muf}-ipdf_{ipdf}")
            
            lo = LO(13600, mur=mur, muf=muf, ipdf=ipdf)

            ptcut = 700.0
            ymax  =   1.7
            for pt in np.geomspace(ptcut, 0.5*sqrts, 400):
                print (pt,
                       lo.dsigma_dpt(pt, 0.0, ymax, include_q=True,  include_g=False),
                       lo.dsigma_dpt(pt, 0.0, ymax, include_q=False, include_g=True))
                
    
