reset
set term pdfcairo enhanced color size 10cm,8cm
set out 'analytics.pdf'
set colors classic

#----------------------------------------------------------------------
# analytic results
#
# We work in the small-angle limit and at lowest orderin alphas
#
# In both cases, the cuts amount to selecting momentum fractions
# between some zmin and some zmax.
#
# We define
#   Iab = \int_zmin^zmax dz Pab(z)
#
# If we start with a quark (or gluon), the fraction of secondary Lund
# planes with a gluon as a leading parton is, respectively
#   f(g|q) = Igq/(Iqq+Igq)
#   f(g|g) = Igg/(Iqg+Igg)
#
# In practice, we have
#
#  Iqq = CF \int_zmin^zmax dz (1+z^2)/(1-z)
#      = CF \int_zmin^zmax dz 2/(1-z) - 1 - z
#      = CF [2 log((1-zmin)/(1-zmax)) - (zmax-zmin) - 1/2 (zmax^2-zmin^2)]
#      = CF [2 log((1-zmin)/(1-zmax)) - (zmax-zmin) (1+(zmax+zmin)/2)]
#
#  Igq = CF \int_zmin^zmax dz (1+(1-z)^2)/z
#      = CF \int_zmin^zmax dz 2/z - 2 + z
#      = CF [2 log(zmax/zmin) - 2 (zmax-zmin) + 1/2 (zmax^2-zmin^2)]
#      = CF [2 log(zmax/zmin) - 2 (zmax-zmin) (1-(zmax+zmin)/4)]
#
#  Iqg = 2 Tf \int_zmin^zmax dz [z^2+(1-z)^2]
#      = 2 Tf/3 [(zmax^3-zmin^3) + ((1-zmin)^3-(1-zmax)^3)]
#      = 2 Tf/3 (zmax-zmin) [2 (zmax^2+zmin zmax+zmin^2) - 3 (zmin+zmax) + 3]
#
#  Igg = 2 CA \int_zmin^zmax dz z/(1-z) + (1-z)/z + z (1-z)
#      = 2 CA \int_zmin^zmax dz 1/(1-z) + 1/z + z (1-z) - 2
#      = 2 CA [log(zmax(1-zmin)/(zmin(1-zmax))) - 2 (zmax-zmin) + (zmax^2-zmin^2)/2 - (zmax^3-zmin^3)/3]

# Quick normalisation check:
# zmax=1/2
#  Iqg = 2 Tf/3 [(1/8-zmin^3) + ((1-zmin)^3-1/8)]
#  Igg = 2 CA [log((1-zmin)/zmin) - 2 (1/2-zmin) + (1/4-zmin^2)/2 - (1/8-zmin^3)/3]
# subtracting the log and zmin->0
#  Iqg = 2 Tf/3 
#  Igg = 2 CA [-11/12]
# sum/(2CA) = -11/12 + TF/(3CA) = (-11CA + 4 Tf)/(12CA)   OK

CF=4.0/3.0
CA=3.0
nf=5.0
TR=0.5
Tf=nf*TR

zmin_def = 0.1
zmax_def = 0.2

Iqq(zmin, zmax) = CF*(2*log((1-zmin)/(1-zmax)) - (zmax-zmin)*(1+(zmax+zmin)/2))
Igq(zmin, zmax) = CF*(2*log(zmax/zmin) - 2*(zmax-zmin)*(1-(zmax+zmin)/4))

Iqg(zmin,zmax) = 2*Tf/3.0 * ((zmax**3-zmin**3) + ((1-zmin)**3-(1-zmax)**3))
Igg(zmin,zmax) = 2*CA*(log(zmax*(1-zmin)/(zmin*(1-zmax))) - 2*(zmax-zmin) + (zmax**2-zmin**2)/2 - (zmax**3-zmin**3)/3)

fg_q(zmin, zmax) = zmin<zmax ? Igq(zmin, zmax)/(Iqq(zmin, zmax)+Igq(zmin, zmax)) : 1/0
fg_g(zmin, zmax) = zmin<zmax ? Igg(zmin, zmax)/(Iqg(zmin, zmax)+Igg(zmin, zmax)) : 1/0



#----------------------------------------------------------------------
# actual plot
set ylabel 'fraction of gluons'
set grid

set title 'Fixing z_{min}=0.1'
set xlabel 'z_{max}'
set xrange [zmin_def:0.5]

plot fg_q(zmin_def, x) w l lt 1 lc 1 lw 2 t 'from quarks', \
     fg_g(zmin_def, x) w l lt 1 lc 3 lw 2 t 'from gluons'

set title 'Fixing z_{max}=0.2'
set xlabel 'z_{min}'
set xrange [0.02:zmax_def]

plot fg_q(x, zmax_def) w l lt 1 lc 1 lw 2 t 'from quarks', \
     fg_g(x, zmax_def) w l lt 1 lc 3 lw 2 t 'from gluons'
     
set title 'both floating'
set xlabel 'z_{min}'
set xrange [0.02:0.5]
set ylabel 'z_{max}'
set yrange [0.02:0.5]
set zlabel 'g fraction'
set zrange [0.6:1]
set contour base
splot fg_q(x, y) w l lt 1 lc 1 lw 2 t 'from quarks', \
      fg_g(x, y) w l lt 1 lc 3 lw 2 t 'from gluons'

# relative difference
set zlabel '[f(g|q)-f(g|g)]/[f(g|q)+f(g|g)]'
set zrange [-0.2:0.05]
splot (fg_q(x, y)-fg_g(x, y))/(fg_q(x, y)+fg_g(x, y)) w l lt 1 lc 1 lw 2 t 'DGLAP'



set out
