reset
set term pdfcairo enhanced color size 10cm,8cm
set out 'Pratios.pdf'
set colors classic

# check that  Pqx/Pgx  increases between z=0 and z=1/2

CF=4.0/3.0
CA=3.0
nf=5.0
TR=0.5
Tf=nf*TR

Pqq(z) = CF*(1+z**2)/(1-z)
Pgq(z) = CF*(1+(1-z)**2)/z
Pqg(z) = 2*Tf*(z**2+(1-z)**2)
Pgg(z) = 2*CA*(z/(1-z) + (1-z)/z + z*(1-z))

#----------------------------------------------------------------------
# actual plot
set ylabel 'P_{qx}/P_{gx}(z)'
set grid

set xlabel 'z'
set xrange [0:0.5]
plot Pqq(x)/Pgq(x) w l lt 1 lc 1 lw 2 t 'x=q', \
     Pqg(x)/Pgg(x) w l lt 1 lc 3 lw 2 t 'x=g'

set out
