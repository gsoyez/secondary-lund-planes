#! /usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0,'../../../PanScales/scripts/')
from hfile import *

pdf_uncert_hess = True
npdf_hess =  40
npdf_mc   = 100

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



Pqq = lambda z: CF*(1+z*z)/(1-z)
Pgq = lambda z: CF*(1+(1-z)*(1-z))/z
Pqg = lambda z: 2*Tf * (np.square(z) + np.square(1-z))
Pgg = lambda z: 2*CA * ((1-z)/z + z/(1-z) + z*(1-z))

Iqq = lambda zmin, zmax: CF*(2*np.log((1-zmin)/(1-zmax)) - (zmax-zmin)*(1+(zmax+zmin)/2))
Igq = lambda zmin, zmax: CF*(2*np.log(zmax/zmin) - 2*(zmax-zmin)*(1-(zmax+zmin)/4))
Iqg = lambda zmin, zmax: 2*Tf/3.0 * ((zmax**3-zmin**3) + ((1-zmin)**3-(1-zmax)**3))
Igg = lambda zmin, zmax: 2*CA*(np.log(zmax*(1-zmin)/(zmin*(1-zmax))) - 2*(zmax-zmin) + (zmax**2-zmin**2)/2 - (zmax**3-zmin**3)/3)

# this assumes zmin<zmax!!
ffg_q = lambda zmin, zmax: Igq(zmin, zmax)/(Iqq(zmin, zmax)+Igq(zmin, zmax)) if zmin<zmax else Pgq(zmin)/(Pqq(zmin)+Pgq(zmin))
ffg_g = lambda zmin, zmax: Igg(zmin, zmax)/(Iqg(zmin, zmax)+Igg(zmin, zmax)) if zmin<zmax else Pgg(zmin)/(Pqg(zmin)+Pgg(zmin))

fg_q = np.vectorize(ffg_q)
fg_g = np.vectorize(ffg_g)


#--------------------------------------------------------------------------------
# processing of NLOJet++ output
#--------------------------------------------------------------------------------

def get_fraction(selection, variable, scale="mur_1-muf_1", ipdf=0):
    ipdf_str = "" if ipdf==0 else f"_{ipdf}"
    tag = "40" if pdf_uncert_hess else "mc"
    dtot = get_array(f"output/3j-PDF4LHC21_{tag}{ipdf_str}-born.ohist", f"sigma_{selection}_v_{variable}-{scale}")
    dg   = get_array(f"output/3j-PDF4LHC21_{tag}{ipdf_str}-born.ohist", f"sigma_{selection}_g_v_{variable}-{scale}")

    f = np.zeros((dtot.shape[0],2))
    f[:,0] = dtot[:,2] if dtot[0,2]>0 else -dtot[:,2]
    ctot = dtot[0,3]
    cg   = dg[0,3]
    f[0,1] = cg/ctot
    for i in range(1,dtot.shape[0]):
        ctot += dtot[i,3]
        cg   += dg  [i,3]
        f[i,1] = cg/ctot

    return f

def get_fraction_and_envelope(selection, variable):
    fcentral = get_fraction(selection, variable)
    f = np.zeros((fcentral.shape[0], 8))
    f[:,0:2] = fcentral

    # scale variation
    f[:,4] = fcentral[:,1]
    f[:,5] = fcentral[:,1]
    for scale in ["mur_0.5-muf_0.5", "mur_0.5-muf_1", "mur_1-muf_0.5", "mur_1-muf_2", "mur_2-muf_1", "mur_2-muf_2"]:
        fvar = get_fraction(selection, variable, scale)
        f[:,4] = np.minimum(f[:,4], fvar[:,1])
        f[:,5] = np.maximum(f[:,5], fvar[:,1])

    # PDF variation
    f[:,6] = 0.0*fcentral[:,1]
    npsd = npdf_hess if pdf_uncert_hess else npdf_mc
    for ipdf in range(1,40):
        fvar = get_fraction(selection, variable, scale)
        f[:,6] += np.square(fvar[:,1]-f[:,1])
    f[:,6] = np.sqrt(f[:,6])
    f[:,7] = f[:,1]-f[:,6]
    f[:,6] = f[:,1]+f[:,6]

    # full variation
    f[:,2] = f[:,1] - np.sqrt(np.square(f[:,4]-f[:,1])+np.square(f[:,6]-f[:,1]))
    f[:,3] = f[:,1] + np.sqrt(np.square(f[:,5]-f[:,1])+np.square(f[:,7]-f[:,1]))
        
    return f

#--------------------------------------------------------------------------------
# plots
#--------------------------------------------------------------------------------

import matplotlib
import matplotlib.pyplot  as plt
import matplotlib.patches as mpatches
from   matplotlib.lines import Line2D as mlines
from   matplotlib.backends.backend_pdf import PdfPages

centrals={"dijet_ptmin" : 150,
          "dijet_ptmax" : 200,
          "dijet_drmin" : 1.0,
          "groom_zmin"  : 0.1,
          "groom_zmax"  : 0.2,
          "groom_drmin" : 0.8,
          }

          

def errorline(ax, legend_patches, f, **kwargs):
    ax.fill_between(f[:,0], f[:,2], f[:,3], lw=0.0, **kwargs, edgecolor=None, alpha=0.2)
    ax.fill_between(f[:,0], f[:,4], f[:,5], lw=0.0, **kwargs, edgecolor=None, alpha=0.2)
    ax.plot        (f[:,0], f[:,1],         lw=1.0, **kwargs)
    rect = mpatches.Patch(facecolor=kwargs['color'], edgecolor=None, alpha=0.2)
    line = mlines([0], [0], **kwargs)
    legend_patches.append((rect,line))

colours=['red', 'green', 'blue', 'magenta']

with PdfPages(f'born-fractions.pdf') as pdf:

    for selection in ["dijet", "groom"]:
        all_vars = ["ptmin", "ptmax", "drmin"] if selection=="dijet" else ["zmin", "zmax", "drmin"]
        all_labs = [r"$p_{t,{{\rm min}}}$", r"$p_{t,{{\rm max}}}$", r"$\Delta R_{{\rm min}}$"] if selection=="dijet" else [r"$z_{{\rm min}}$", r"$z_{{\rm max}}$", r"$\Delta R_{{\rm min}}$"]

        for v,l in zip(all_vars, all_labs):
            fig, ax = plt.subplots(figsize=(4.0, 3))

            ax.set_title(f'{selection} selection - '+l+' dependence')

            ax.grid(True, which='both', lw=0.5, ls=':', zorder=0)
            ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
        
            ax.set_xlabel(l)
            if selection=="dijet":
                # zmin=ptmin/(ptlead+ptmin)
                # zmax=ptmax/(ptlead+ptmax)
                if v=="ptmin":
                    ax.set_xlim(100.0, 200.0)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = fg_q(an_x/(700.0+an_x), 200.0/(200.0+700.0))
                    an_g = fg_g(an_x/(700.0+an_x), 200.0/(200.0+700.0))
                elif v=="ptmax":
                    ax.set_xlim(150.0, 300.0)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = fg_q(100.0/(700.0+100.0), an_x/(an_x+700.0))
                    an_g = fg_g(100.0/(700.0+100.0), an_x/(an_x+700.0))
                else:
                    ax.set_xlim(0.4, 1.2)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = 0.0*an_x+fg_q(100.0/(700.0+100.0), 200.0/(200.0+700.0))
                    an_g = 0.0*an_x+fg_g(100.0/(700.0+100.0), 200.0/(200.0+700.0))
            else:
                if v=="zmin":
                    ax.set_xlim(0.05, 0.2)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = fg_q(an_x,0.2)
                    an_g = fg_g(an_x,0.2)
                elif v=="zmax":
                    ax.set_xlim(0.1, 0.5)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = fg_q(0.1,an_x)
                    an_g = fg_g(0.1,an_x)
                else:
                    ax.set_xlim(0.0, 1.2)
                    an_x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 51)
                    an_q = 0.0*an_x+fg_q(0.1,0.2)
                    an_g = 0.0*an_x+fg_g(0.1,0.2)
            #ax.set_xlim(0.0, 5.0)
            #xt=np.linspace(0,5,6)
            #ax.set_xticks(xt)
            #ax.set_xticklabels([f'{x}' for x in xt])

            ax.set_ylabel("gluon fraction")

            legend_patches, legend_labels = ax.get_legend_handles_labels()
            
            ax.set_ylim(0.83,0.93)
            yt=np.linspace(0.83,0.93,11)
            ax.set_yticks(yt)
            ax.set_yticklabels([f'{x:g}' for x in yt])

            # NLOJet
            f = get_fraction_and_envelope(selection, v)
            errorline(ax, legend_patches, f, ls='-', color='blue')
            legend_labels.append(r"Exact $\mathcal{O}(\alpha_s)$")
            
            # analytics
            ax.plot(an_x, an_q, ls='-', color='red',   label='collinear quark')
            ax.plot(an_x, an_g, ls='-', color='green', label='collinear gluon')

            y1=0.125
            y2=0.045
            if selection=="dijet":
                ax.text(0.05, y1, r'anti-$k_t(R=0.4)$, $p_{t,{\rm lead}}>700$ GeV, $|y|<1.7$',  ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

                if v=="ptmin":
                    ax.text(0.05, y2, r'$p_{t,{\rm sublead}}<200$ GeV, $1<\Delta R<1.2$',    ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                elif v=="ptmax":
                    ax.text(0.05, y2, r'$150<p_{t,{\rm sublead}}$ GeV, $1<\Delta R<1.2$',   ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                else:
                    ax.text(0.05, y2, r'$150<p_{t,{\rm sublead}}<200$ GeV, $\Delta R<1.2$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
            else:
                ax.text(0.05, y1, r'anti-$k_t(R=1.2)$, $p_{t,{\rm lead}}>1$ TeV, $|y|<1.7$',   ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

                if v=="zmin":
                    ax.text(0.05, y2, r'$z_g<0.2$, $\Delta R|y|>0.8$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                elif v=="zmax":
                    ax.text(0.05, y2, r'$0.1<z_g$, $\Delta R|y|>0.8$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                else:
                    ax.text(0.05, y2, r'$0.1<z_g<0.2$',                ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

            p,l = ax.get_legend_handles_labels()
            for pp in p: legend_patches.append(pp)
            for ll in l: legend_labels.append(ll)
            ax.legend(legend_patches, legend_labels)
            pdf.savefig(bbox_inches='tight')
            plt.close()



# 
# reset
# set term pdfcairo enhanced color
# set out 'born-fractions.pdf'
# 
# set ylabel 'gluon fraction'
# set grid
# 
# set title 'dijet selection - p_{t,max} dependence'
# set xlabel 'p_{t,max} [GeV]'
# 
# plot '< mergeidx.pl -f 
# 
# set out
# 
