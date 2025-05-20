#! /usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0,'../../../PanScales/scripts/')
from hfile import *

pdf_uncert_hess = False
npdf_hess =  40
npdf_mc   = 100

#--------------------------------------------------------------------------------
# processing of NLOJet++ output
#--------------------------------------------------------------------------------

def get_fraction(selection, variable, scale="mur_1-muf_1", ipdf=0):
    ipdf_str = "" if ipdf==0 else f"_{ipdf}"
    tag = "40" if pdf_uncert_hess else "mc"
    dtot = get_array(f"output/3j-smallr-PDF4LHC21_{tag}{ipdf_str}-born.ohist", f"sigma_{selection}_v_{variable}-{scale}")
    dg   = get_array(f"output/3j-smallr-PDF4LHC21_{tag}{ipdf_str}-born.ohist", f"sigma_{selection}_g_v_{variable}-{scale}")

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
    # npsd = npdf_hess if pdf_uncert_hess else npdf_mc
    # for ipdf in range(1,40):
    #     fvar = get_fraction(selection, variable, scale)
    #     f[:,6] += np.square(fvar[:,1]-f[:,1])
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

ptcut=700.0
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

with PdfPages(f'born-fractions-smallr.pdf') as pdf:

    for selection in ["dijet", "groom"]:
        all_vars = ["ptmin", "ptmax", "drmin"] if selection=="dijet" else ["zmin", "zmax", "drmin"]
        all_labs = [r"$p_{t,{{\rm min}}}$ [GeV]", r"$p_{t,{{\rm max}}}$ [GeV]", r"$\Delta R_{{\rm min}}$"] if selection=="dijet" else [r"$z_{{\rm min}}$", r"$z_{{\rm max}}$", r"$\Delta R_{{\rm min}}$"]

        for v,l in zip(all_vars, all_labs):
            print (f"{selection=}, {v=}")
            fig, ax = plt.subplots(figsize=(4.0, 3))

            ax.set_title(f'{selection} selection - '+l.removesuffix(" [GeV]")+' dependence')

            ax.grid(True, which='both', lw=0.5, ls=':', zorder=0)
            ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
        
            ax.set_xlabel(l)

            analytics = get_array("../analytics/analytic-born-fractions-smallr.res",
                                  f"sel_{selection}-dep_{v}")
            an_x = analytics[:,0]
            an_q = analytics[:,1]
            an_g = analytics[:,2]
            if selection=="dijet":
                if v=="ptmin":   ax.set_xlim(100.0, 200.0)
                elif v=="ptmax": ax.set_xlim(150.0, 300.0)
                else:            ax.set_xlim(0.2, 0.6)
            
            else:
                if v=="zmin":    ax.set_xlim(0.05, 0.2)
                elif v=="zmax":  ax.set_xlim(0.1, 0.5)
                else:            ax.set_xlim(0.0, 0.4)
                
            ax.set_ylabel("gluon fraction")

            legend_patches, legend_labels = ax.get_legend_handles_labels()
            
            ax.set_ylim(0.798,0.962)
            yt=np.linspace(0.8,0.96,9)
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
                ax.text(0.035, y1, r'anti-$k_t(R=0.2)$, $p_{t,{\rm lead}}>700$ GeV, $|y|<1.7$',  ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

                if v=="ptmin":
                    ax.text(0.035, y2, r'$p_{t,{\rm sublead}}<200$ GeV, $0.4<\Delta R<0.6$',    ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                elif v=="ptmax":
                    ax.text(0.035, y2, r'$150<p_{t,{\rm sublead}}$ GeV, $0.4<\Delta R<0.6$',   ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                else:
                    ax.text(0.035, y2, r'$150<p_{t,{\rm sublead}}<200$ GeV, $\Delta R<0.6$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
            else:
                ax.text(0.035, y1, r'anti-$k_t(R=0.4)$, $p_{t,{\rm lead}}>1$ TeV, $|y|<1.7$',   ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

                if v=="zmin":
                    ax.text(0.035, y2, r'$z_g<0.2$, $\Delta R>0.3$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                elif v=="zmax":
                    ax.text(0.035, y2, r'$0.1<z_g$, $\Delta R>0.3$', ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')
                else:
                    ax.text(0.035, y2, r'$0.1<z_g<0.2$',             ha='left', va='baseline', transform=ax.transAxes, fontsize=9, color='black')

            p,l = ax.get_legend_handles_labels()
            for pp in p: legend_patches.append(pp)
            for ll in l: legend_labels.append(ll)
            ax.legend(legend_patches, legend_labels, loc='upper right')
            pdf.savefig(bbox_inches='tight')
            plt.close()


