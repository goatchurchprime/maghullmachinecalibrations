import sys
sys.path.append("/home/goatchurch/geom3d/maghullmachinecalibrations")
from basicgeo import P2, I1
from math import sqrt, acos, sin, cos, degrees, radians
import scipy.optimize

Xfake = (744, 746, 747, 694.5, 693.0, 584.6, 0.0707, 518.4, 387.9)
j0range = I1(-300, 100)
j1range = I1(-230, 95)

def abanglefromtriangleabc(a, b, c):
    ac = (a*a + b*b - c*c)/(2*a*b)
    #assert -1<=ac<=1, ac
    return acos(max(-1,min(1,ac))) 

def jointstocartesian(j0, j1, Xmodel):
    a, b, c, d, e, f, ah, tblx, armx = Xmodel
    ab = abanglefromtriangleabc(a, b, c) 
    ac = abanglefromtriangleabc(a, c, b) 
    tbl = j0 + tblx
    arm = j1 + armx

    aTBL = abanglefromtriangleabc(b, d, tbl)
    aARM = abanglefromtriangleabc(c, e, arm) 

    avec = ab - aTBL 
    pvec = avec + ac - aARM + ah 

    return P2(0 - a*sin(avec) + f*sin(pvec),  
              d - a*cos(avec) + f*cos(pvec))


Fa, Fb, Fc, Fd, Fe, Ff, Fah, Ftblx, Farmx = Xfake

# proof of calculation from 4 dowel holes with a fixed j0
j0 = j0range.Along(0.5)
j1s = [j1range.Along(lam)  for lam in [0.1, 0.4, 0.6, 0.9]]
DHs = [ jointstocartesian(j0, j1, Xfake)  for j1 in j1s ]

# these are the measurements between the 4 dowel holes
d01 = (DHs[0] - DHs[1]).Len()
d02 = (DHs[0] - DHs[2]).Len()
d12 = (DHs[1] - DHs[2]).Len()
d03 = (DHs[0] - DHs[3]).Len()

def circumcircleradius(a, b, c):
    s = (a+b+c)/2
    areasq = s*(s-a)*(s-b)*(s-c)
    return a*b*c/(4*sqrt(areasq))
Cf = circumcircleradius(d01, d02, d12)
print("f:", Cf, "from", Ff)

sweepang01 = acos(1-d01**2/(2*Cf*Cf))
sweepang02 = acos(1-d02**2/(2*Cf*Cf))
sweepang03 = acos(1-d03**2/(2*Cf*Cf))

# these sweep angles rely on the three unknowns e, c, armx
# so we can solve by minimizing across three equations (dimensions)
def fun1(X):
    e, c, armx = X
    absangsX = [ acos((e**2 + c**2 - (armx + j1)**2)/(2*e*c))  for j1 in j1s ]
    sa01X = absangsX[1] - absangsX[0]
    sa02X = absangsX[2] - absangsX[0]
    sa03X = absangsX[3] - absangsX[0]
    return (sa01X-sweepang01)**2 + (sa02X-sweepang02)**2 + (sa03X-sweepang03)**2

# this minimization converges from almost any starting point
res = scipy.optimize.minimize(fun=fun1, x0=(Cf, Cf, Cf), method='Nelder-Mead', options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
Ce, Cc, Carmx = res.x
print("e:", Ce, "from", Fe)
print("c:", Cc, "from", Fc)
print("armx:", Carmx, "from", Farmx)

# the third side of triangle af is the radius for fixed j1
# this depends on the new values a, b, ah
# Needing a, b, ah in addition to already found Cc, Ce, Carmx, Cf

# this will require at least three radii to be measured from 9 dowels
#arm = j1 + armx
#ac = abanglefromtriangleabc(a, c, b)
#aARM = abanglefromtriangleabc(c, e, arm)
#af = ac - aARM + ah
#Rf = sqrt(a**2 + f**2 - 2*a*f*cos(af))

# measure the three radii by swinging the table in j0 for the first three j1s
j0s = [j0range.Along(lam)  for lam in [0.2, 0.5, 0.8]]
Radforj1s = [ ]
for j1 in j1s:
    EHs = [ jointstocartesian(j0, j1, Xfake)  for j0 in j0s ]
    e01 = (EHs[0] - EHs[1]).Len()
    e02 = (EHs[0] - EHs[2]).Len()
    e12 = (EHs[1] - EHs[2]).Len()
    Rf = circumcircleradius(e01, e02, e12)
    Radforj1s.append((j1, Rf))

# values Cc, Ce, Carmx, Cf already known
# solved are a and acah
Facah = abanglefromtriangleabc(Fa, Fc, Fb) + Fah
def Radforj1(j1, a, acah):
    arm = j1 + Carmx
    aARM = abanglefromtriangleabc(Cc, Ce, arm)
    af = acah - aARM
    Rf = sqrt(a**2 + Cf**2 - 2*a*Cf*cos(af))
    return Rf

def fun2(X):
    a, acah = X
    return sum((Radforj1(j1, a, acah) - Rf)**2  for j1, Rf in Radforj1s[:3])

# another good convergence due to small number of unknowns solved
# but needs at least 3 radii of table swings to work
res = scipy.optimize.minimize(fun=fun2, x0=(Cc, 0.1), method='Nelder-Mead', options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
Ca, Cacah = res.x
print("a:", Ca, "from", Fa)
print("acah:", Cacah, "from", Facah)

# remaining unknowns are d, c, ah, tblx (minus acah)

# consider the triangle d, b, tbl and 4 dowels at a fixed j1 as before
j1 = j1range.Along(0.5)
j0s = [ j1range.Along(lam)  for lam in [0.1, 0.4, 0.6, 0.9] ]
EHs = [ jointstocartesian(j0, j1, Xfake)  for j0 in j0s ]
e01 = (EHs[0] - EHs[1]).Len()
e02 = (EHs[0] - EHs[2]).Len()
e03 = (EHs[0] - EHs[3]).Len()

Rf = Radforj1(j1, Ca, Cacah)
Tsweepang01 = acos(1-e01**2/(2*Rf*Rf))
Tsweepang02 = acos(1-e02**2/(2*Rf*Rf))
Tsweepang03 = acos(1-e03**2/(2*Rf*Rf))

def fun3(X):
    d, b, tblx = X
    absangsX = [ acos((d**2 + b**2 - (tblx + j0)**2)/(2*d*b))  for j0 in j0s ]
    sa01X = absangsX[1] - absangsX[0]
    sa02X = absangsX[2] - absangsX[0]
    sa03X = absangsX[3] - absangsX[0]
    return (sa01X-Tsweepang01)**2 + (sa02X-Tsweepang02)**2 + (sa03X-Tsweepang03)**2

# another minimization on this triangle should have same convergence stability
res = scipy.optimize.minimize(fun=fun3, x0=(Cc, Cc, Farmx), method='Nelder-Mead', options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
Cd, Cb, Ctblx = res.x
print("d:", Cd, "from", Fd)
print("c:", Cc, "from", Fc)
print("tblx:", Ctblx, "from", Ftblx)

# finally get ah back out
Cah = Cacah - abanglefromtriangleabc(Ca, Cc, Cb)
print("ah:", Cah, "from", Fah)

# so that's 11 dowel holes in this configuration
#   
#   00     01     02
#   
#   
#   10     11     12     13
#
#
#   20     21     22
#
#
#          31
#
# where all horizontal distances in the 3x3 square are required (18 measurements)
# in addition to 01-31 and 10-13, giving a total of 20
