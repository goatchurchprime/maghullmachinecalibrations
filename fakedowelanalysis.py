import sys
sys.path.append("/home/goatchurch/geom3d/maghullmachinecalibrations")
from doesmachinemodel import jointstocartesian, cartesiantojoints
from doesmachinemodel import Xsample, j0range, j1range
from basicgeo import P2, I1
from math import sqrt, acos, cos, degrees, radians

# generate the joint positions for the dowel holes
j0s = [j0range.Along(lam)  for lam in [0.2, 0.5, 0.8]]
j1s = [j1range.Along(lam)  for lam in [0.2, 0.5, 0.8]]
jarray = [[(j0, j1)  for j1 in j1s]  for j0 in j0s]

# select the joint pairs where distances are measured
jointpairs = [ ]
for j0 in j0s:
    for j1lo in j1s:
        for j1hi in j1s:
            if j1lo < j1hi:
                jointpairs.append(((j0,j1lo),(j0,j1hi)))
for j1 in j1s:
    for j0lo in j0s:
        for j0hi in j0s:
            if j0lo < j0hi:
                jointpairs.append(((j0lo,j1),(j0hi,j1)))

sendactivity("clearallpoints")
sendactivity("points", points=[j[0]  for j in jointpairs])

Xfake = Xsample
# convert to measurement distances
dm = { }
for ((j0,j1),(k0,k1)) in jointpairs:
    dm[((j0,j1),(k0,k1))] = (jointstocartesian(j0, j1, Xfake) - jointstocartesian(k0, k1, Xfake)).Len()

def circumcircleradius(a, b, c):
    s = (a+b+c)/2
    areasq = s*(s-a)*(s-b)*(s-c)
    return a*b*c/(4*sqrt(areasq))
    
j0 = j0s[0]
d01 = dm[((j0,j1s[0]),(j0,j1s[1]))]    
d12 = dm[((j0,j1s[1]),(j0,j1s[2]))]    
d02 = dm[((j0,j1s[0]),(j0,j1s[2]))]    

Cf = circumcircleradius(d01,d12,d02)
print(Cf)

ang01 = acos(1-d01**2/(2*Cf*Cf))
ang12 = acos(1-d12**2/(2*Cf*Cf))
ang02 = ang01 + ang12

a, b, c, d, E, f, ah, tblx, armx = Xfake
angc0 = (acos((E**2 + c**2 - (armx + j1s[0])**2)/(2*E*c)))
angc1 = (acos((E**2 + c**2 - (armx + j1s[1])**2)/(2*E*c)))
angc2 = (acos((E**2 + c**2 - (armx + j1s[2])**2)/(2*E*c)))
a01 = angc1 - angc0
a12 = angc2 - angc1
angc2 - angc0

# ang01, ang02, i1, i2 are known, but e, c, j are not

import scipy.optimize

i1 = j1s[1] - j1s[0]
i2 = j1s[2] - j1s[0]
def fun1(X):
    E, c, j = X
    ce2 = 2*E*c
    Da01 = acos((E**2 + c**2 - (j+i1)**2)/ce2) - acos((E**2 + c**2 - j**2)/ce2)
    Da02 = acos((E**2 + c**2 - (j+i2)**2)/ce2) - acos((E**2 + c**2 - j**2)/ce2)
    return (Da01 - ang01)**2 + (Da02 - ang02)**2
fun1((E+1, c, armx + j1s[0]))
j = armx + j1s[0]
(E, c, armx + j1s[0])

scipy.optimize.minimize(fun=fun1, x0=(E, c+1, armx + j1s[0]), method='Nelder-Mead', options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})

sendactivity("clearalltriangles")

# explore the linearity of the table drive
acc = [ ]
acs = [ ]
for k in range(1, 1000):
    i = I1(-armx+c-E, -armx+c+E).Along(k/1000)
    Da = acos((E**2 + c**2 - (armx+i)**2)/ce2) - acos((E**2 + c**2 - armx**2)/ce2)
    acc.append((i, degrees(Da)))
    Das = 2*(armx+i)/sqrt(1-((E**2 + c**2 - (armx+i)**2)/ce2)**2)/ce2
    acs.append((i, degrees(Das)))
sendactivity("contours", contours=[acc])
sendactivity("contours", contours=[acs])
sendactivity("contours", contours=[[acc[0], acc[-1]]], materialnumber=1)
                                  



