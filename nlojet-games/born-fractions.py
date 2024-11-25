#! /usr/bin/env python3

import numpy as np
import sys
sys.path.insert(0,'../../../PanScales/scripts/')
from hfile import *

pdf_uncert_hess = True
npdf_hess =  40
npdf_mc   = 100

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

with PdfPages(f'born-fractions.pdf') as pdf:

    for selection in ["dijet", "groom"]:
        all_vars = ["ptmin", "ptmax", "drmin"] if selection=="dijet" else ["zmin", "zmax", "drmin"]
        all_labs = [r"$p_{t,{{\rm min}}}$", r"$p_{t,{{\rm max}}}$", r"$\Delta R_{{\rm min}}$"] if selection=="dijet" else [r"$z_{{\rm min}}$", r"$z_{{\rm max}}$", r"$\Delta R_{{\rm min}}$"]

        for v,l in zip(all_vars, all_labs):
            print (f"{selection=}, {v=}")
            fig, ax = plt.subplots(figsize=(4.0, 3))

            ax.set_title(f'{selection} selection - '+l+' dependence')

            ax.grid(True, which='both', lw=0.5, ls=':', zorder=0)
            ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
        
            ax.set_xlabel(l)

            analytics = get_array("../analytics/analytic-born-fractions.res",
                                  f"sel_{selection}-dep_{v}")
            an_x = analytics[:,0]
            an_q = analytics[:,1]
            an_g = analytics[:,2]
            an0_q = None
            an0_g = None
            do_an_band=False
            if selection=="dijet":
                if v=="ptmin":   ax.set_xlim(100.0, 200.0)
                elif v=="ptmax": ax.set_xlim(150.0, 300.0)
                else:            ax.set_xlim(0.4, 1.2)
                # uncomment the next 2 lines to get the fraction under
                # the assumotion that the leading jet has exactly
                # pt=ptmin                
                #an0_q = analytics[:,-2]
                #an0_g = analytics[:,-1]

                # uncomment the next lines to get a band
                # do_an_band=True
                # an_q = np.zeros((an_x.shape[0],6))
                # an_g = np.zeros((an_x.shape[0],6))
                # an_q[:,0] = analytics[:,0]
                # an_q[:,1] = analytics[:,1]
                # an_q[:,2] = analytics[:,1]+analytics[:,3]
                # an_q[:,3] = analytics[:,1]+analytics[:,5]
                # an_q[:,4] = analytics[:,1]+analytics[:,7]
                # an_q[:,5] = analytics[:,1]+analytics[:,9]
                # an_g[:,0] = analytics[:,0]
                # an_g[:,1] = analytics[:,2]
                # an_g[:,2] = analytics[:,2]+analytics[:,4]
                # an_g[:,3] = analytics[:,2]+analytics[:,6]
                # an_g[:,4] = analytics[:,2]+analytics[:,8]
                # an_g[:,5] = analytics[:,2]+analytics[:,10]
            
            else:
                if v=="zmin":    ax.set_xlim(0.05, 0.2)
                elif v=="zmax":  ax.set_xlim(0.1, 0.5)
                else:            ax.set_xlim(0.0, 1.2)
                
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
            if do_an_band:
                errorline(ax, legend_patches, an_q, ls='-', color='red')
                legend_labels.append(r"collinear quark")
                errorline(ax, legend_patches, an_g, ls='-', color='green')
                legend_labels.append(r"collinear gluon")
            else:
                ax.plot(an_x, an_q, ls='-', color='red',   label='collinear quark')
                ax.plot(an_x, an_g, ls='-', color='green', label='collinear gluon')
                if an0_q is not None:
                    ax.plot(an_x, an0_q, ls='--', color='red')
                    ax.plot(an_x, an0_g, ls='--', color='green')

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
