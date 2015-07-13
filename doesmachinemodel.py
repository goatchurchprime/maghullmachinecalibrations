from math import sqrt, atan, cos, sin, acos, radians, atan2
from basicgeo import P2, I1

Xsample = (744, 746, 747, 694.5, 693.0, 584.6, 0.0707, 518.4, 387.9)
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

def clengthfromtriangleabangle(a, b, ang):
    return sqrt(a*a + b*b - 2*a*b*cos(ang)) 
def cartesiantojoints(p, Xmodel):
    a, b, c, d, e, f, ah, tblx, armx = Xmodel
    ab = abanglefromtriangleabc(a, b, c) 
    ac = abanglefromtriangleabc(a, c, b) 
    x, y = p
    gx = x
    gy = y - d
    g = sqrt(gx*gx + gy*gy)

    avec = abanglefromtriangleabc(a, g, f) - atan2(gx, -gy) 
    af = abanglefromtriangleabc(a, f, g) 
    pvec = avec + af 
    aTBL = ab - avec 
    aARM = ac - af + ah

    tbl = clengthfromtriangleabangle(b, d, aTBL) 
    arm = clengthfromtriangleabangle(c, e, aARM) 

    return (tbl - tblx, arm - armx)

