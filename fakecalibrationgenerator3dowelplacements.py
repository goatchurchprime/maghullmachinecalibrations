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


j0range, j1range = (-240, 140), (-210, 80)

# fake data generator
Xf = (744, 746, 747, 694.5, 693.0, 584.6, 0.0707, 518.4, 387.9)  # fake

n0, n1 = 3, 3
def Along(lam, lo, hi):  return lo*(1-lam) + hi*lam
j01 = [ (Along(i0/n0, j0range[0], j0range[1]), Along(i1/n1, j1range[0], j1range[1]))  for i0 in range(n0+1)  for i1 in range(n1+1) ]
jpairs = [ (j01[i], j01[i1])  for i in range(len(j01))  for i1 in range(i+1, len(j01)) ]

def jprdists(X):
    def djpr(jA, jB):
        pA, pB = procforward(jA[0], jA[1], X), procforward(jB[0], jB[1], X)
        return sqrt((pA[0]-pB[0])**2 + (pA[1]-pB[1])**2)
    return [ djpr(jA, jB)  for jA, jB in jpairs ]
Jdistf = jprdists(Xf)

def fun1(Xs):
    X = insertXfsub(Xs)
    return sum((jd0-jd1)**2  for jd0, jd1 in zip(jprdists(X), Jdistf))


def knockXfsub(X):
    return [X[i]  for i in Xfsub]
def insertXfsub(X):
    return [ i not in Xfsub and Xf[i] or X[Xfsub.index(i)]  for i in range(9) ]

print(fun1(Xf))


method = 'Powell'
method = 'Nelder-Mead'
Xfsub = [0,1,2,3,5,6,8]#[0,4,5,7,8] # list(range(6))
Xfk = knockXfsub(Xf)
X0k = tuple(x-0.05  for x in Xfk)
bnds = [ (x-5,x+5)  for x in X0k ]
tstart = time.time()
res = scipy.optimize.minimize(fun=fun1, x0=X0k, bounds=bnds, method=method, options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
Xoptk = tuple(res.x)
Xopt = insertXfsub(Xoptk)
print("Elapsed seconds %.0f  orig"%(time.time() - tstart), fun1(X0k), "opt", fun1(Xoptk))
print(", ".join("%.3f"%(x-x0)  for x, x0 in zip(Xoptk, Xfk)))

d = 0.0001
def fun1d(Xk, i):
    lXk0, lXk1 = list(Xk[:]), list(Xk[:])
    lXk0[i], lXk1[i] = Xk[i]-d, Xk[i]+d
    return fun1(lXk0)-mo, fun1(lXk1)-mo
Xk = list(res.x)
mo = fun1(Xk)
print([fun1d(Xk, i)  for i in range(len(Xk))])

res = scipy.optimize.minimize(fun=fun1, x0=Xk, bounds=bnds, method="BFGS", options={"xtol":0.00000001, "ftol":1e-20, "maxiter":9000})
print(res.fun)
Xk = list(res.x)

print(fun1(Xk))
Xk[0]-=d

print(fun1(Xoptk))

