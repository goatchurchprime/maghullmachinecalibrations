from math import sqrt, atan, cos, sin, acos, radians, atan2
import os, time, scipy.optimize

try:  sendactivity
except NameError:
    def sendactivity(*args, **kwargs):  pass
    
R2 = 132.525
X0 = (744.57703517328832, 746.19009986079459, 741.95547170764314, 695.42905812855656, 693.07538973655983, 582.64403220610257, 0.0764480526, 690, 400)

def abanglefromtriangleabc(a, b, c):
    ac = (a*a + b*b - c*c)/(2*a*b)
    #assert -1<=ac<=1, ac
    return acos(max(-1,min(1,ac))) 
def procforward(j0, j1, X, diskcen=(0,0)):
    a, b, c, d, e, f, ah, tblx, armx = X
    ab = abanglefromtriangleabc(a, b, c) 
    ac = abanglefromtriangleabc(a, c, b) 
    tbl = j0 + tblx
    arm = j1 + armx

    aTBL = abanglefromtriangleabc(b, d, tbl)
    aARM = abanglefromtriangleabc(c, e, arm) 

    avec = ab - aTBL 
    pvec = avec + ac - aARM + ah 

    dx, dy = diskcen
    return (0 - a*sin(avec) + f*sin(pvec) - dx,  
            d - a*cos(avec) + f*cos(pvec) - dy)

def clengthfromtriangleabangle(a, b, ang):
    return sqrt(a*a + b*b - 2*a*b*cos(ang)) 
def procback(x, y, X):
    a, b, c, d, e, f, ah, tblx, armx = X
    ab = abanglefromtriangleabc(a, b, c) 
    ac = abanglefromtriangleabc(a, c, b) 
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

# fake data generator
Cf = (231.8, 291.1)  # centrepoint
R2 = 132.525
Xf = (744, 746, 747, 694.5, 693.0, 584.6, 0.0707, 518.4, 387.9)  # fake
j1values = [-213.0, -183.1, -153.2, -123.3, -93.4, -63.5, -33.6, -3.7, 26.2, 56.1, 86.0]
def R2circdist(j0, j1):
    x, y = procforward(j0, j1, Xf)
    r = sqrt((Cf[0] - x)**2 + (Cf[1] - y)**2)
    return (r - R2)
def R2circcut(j0lo, j0hi, j1):
    assert R2circdist(j0lo, j1) < 0 < R2circdist(j0hi, j1)
    while abs(j0hi - j0lo) > 1e-10:
        j0mid = (j0lo + j0hi)/2
        if R2circdist(j0mid, j1) < 0:
            j0lo = j0mid
        else:
            j0hi = j0mid
    return j0mid

def MakeFakeDiskSamplesBisect():
    return [(R2circcut(0, 300, j1), j1)  for j1 in j1values] + \
           [(R2circcut(0, -300, j1), j1)  for j1 in j1values]

def BestCentreR2(pts, R2):
    cx0, cy0 = sum(x  for x, y in pts)/len(pts), sum(y  for x, y in pts)/len(pts)
    def centrequality(c):
        dss = [ sqrt((c[0] - x)**2 + (c[1] - y)**2)  for x, y in pts ]
        return sum((d - R2)**2  for d in dss)/len(pts)
    b = 25
    bnds = ((cx0-b,cx0+b), (cy0-b,cy0+b))
    res = scipy.optimize.minimize(fun=centrequality, x0=(cx0, cy0), bounds=bnds)
    return tuple(res.x), centrequality(res.x)

# main code
disksamples = MakeFakeDiskSamplesBisect()

pts = [ procforward(j0, j1, Xf)  for j0, j1 in disksamples ]
c, q = BestCentreR2(pts, R2)

def transd(p, c, d):
    vx, vy = p[0], p[1] - d
    ux, uy = c[0], c[1] - d
    ulen = sqrt(ux**2 + uy**2)
    return (vx*ux + vy*uy)/ulen, (vx*uy - vy*ux)/ulen
def fun1(Xs):
    X = insertXfsub(Xs)
    pts = [ procforward(j0, j1, X)  for j0, j1 in disksamples ]
    c, q = BestCentreR2(pts, R2)
    if pfac == 0:  return q
    ptax = [transd(p, c, X[3])  for p in pts]
    ax = sum((pa0[0]-pa1[0])**2 + (pa0[1]+pa1[1])**2  for pa0, pa1 in zip(ptax[::2], ptax[1::2]))
    return q + pfac*ax

def radrange(X):
    pts = [ procforward(j0, j1, X)  for j0, j1 in disksamples ]
    c, q = BestCentreR2(pts, R2)
    dss = [ sqrt((c[0] - x)**2 + (c[1] - y)**2)  for x, y in pts ]
    ptax = [transd(p, c, X[3])  for p in pts]
    axl = [(abs(pa0[0]-pa1[0]), abs(pa0[1]+pa1[1]))  for pa0, pa1 in zip(ptax[::2], ptax[1::2])]
    print("minr", min(dss), "maxr", max(dss), "dalong", max(ax[0]  for ax in axl), "dside", max(ax[1]  for ax in axl))


def knockXfsub(X):
    return [X[i]  for i in Xfsub]
def insertXfsub(X):
    return [ i not in Xfsub and Xf[i] or X[Xfsub.index(i)]  for i in range(9) ]

method = 'Powell'
Xfsub = [5,6,7,8]#[0,4,5,7,8] # list(range(6))
Xfk = knockXfsub(Xf)
X0k = tuple(x+0.1  for x in Xfk)
bnds = [ (x-5,x+5)  for x in X0k ]
pfac = 0*0.08  # control this factor for how much input made by axis symmetry
tstart = time.time()
res = scipy.optimize.minimize(fun=fun1, x0=X0k, bounds=bnds, method=method, options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
Xoptk = tuple(res.x)
Xopt = insertXfsub(Xoptk)
pts = [ procforward(j0, j1, Xopt)  for j0, j1 in disksamples ]
c, q = BestCentreR2(pts, R2)
print("Elapsed seconds %.0f  orig"%(time.time() - tstart), fun1(X0k), "opt", fun1(Xoptk))
print(", ".join("%.2f"%(x-x0)  for x, x0 in zip(Xoptk, Xfk)), c[0]-Cf[0], c[1]-Cf[1])

