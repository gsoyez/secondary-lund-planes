# gg -> gg
#---------
def Mgg2gg(s, t, u):
    return 0.5 * 9.0/2.0*(3.0-(u*t)/(s*s)-(u*s)/(t*t)-(s*t)/(u*u))

# gg -> qqbar
#------------
def Mgg2qqbar(s, t, u):
    Nf=5
    return Nf*(1.0/6.0*(u/t+t/u) - 3.0/8.0*(u*u+t*t)/(s*s))

# qg -> qg
#---------
def Mqg2qg(s, t, u):
    return (u*u+s*s)/(t*t) - 4.0/9.0*(u/s+s/u)

# qq' -> qq'
#-----------
def Mqqp2qqp(s, t, u):
    return 4.0/9.0*(s*s+u*u)/(t*t)

# qqbar -> q'qbar'
#-----------------
def Mqqbar2qpqbarp(s, t, u):
    Nf=5
    return (Nf-1)*4.0/9.0*(t*t+u*u)/(s*s)

# qq -> qq
#---------
def Mqq2qq(s, t, u):
    return 0.5*(4.0/9.0*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u)) - 8.0/27.0*(s*s)/(t*u))

# qqbar -> qqbar
#---------------
def Mqqbar2qqbar(s, t, u):
    return 4.0/9.0*((t*t+u*u)/(s*s)+(s*s+u*u)/(t*t)) - 8.0/27.0*(u*u)/(s*t)

# qqbar -> gg
#------------
def Mqqbar2gg(s, t, u):
    return 0.5*(32.0/27.0*(u/t+t/u) - 8.0/3.0*(t*t+u*u)/(s*s))
