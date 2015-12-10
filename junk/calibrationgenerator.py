from math import sqrt, atan, cos, sin, acos, radians
import scipy.optimize
import os, time

# circle sizes = 132.525, 79.69, 17.055

try:  sendactivity
except NameError:
    def sendactivity(*args, **kwargs):  pass
    
fname = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/probing2015-05-22a.txt"


X0 = (744.57703517328832, 746.19009986079459, 741.95547170764314, 695.42905812855656, 693.07538973655983, 582.64403220610257, 0.0764480526, 690, 400)
#X0 = (743.3133003594869, 756.82976800937217, 743.43834112653599, 698.16794077213683, 688.61101385919824, 580.7463415098855, 0.097855119579137784, 518.65650220023929, 408.6731700849632)

halsampleHz = 1000
Nmergablegap = 160  # when finding samples in contact
aXoffset = 0#-0.4837

def abanglefromtriangleabc(a, b, c):
    return acos((a*a + b*b - c*c)/(2*a*b)) 
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

class ContactSequence:
    def __init__(self, i, sampleline, nprevgap, N):
        self.i = i
        self.samplelines = [ sampleline ]
        self.nprevgap = nprevgap
        self.ninternalgap = 0
        self.N = N
    def ApplyConversion(self, X):
        self.pts = [ procforward(float(sampleline[0]), float(sampleline[1]), X)  for sampleline in self.samplelines ]
    def guessisel(self, mx, my):
        dseq = [((x - mx)**2 + (y - my)**2, i)  for i, (x, y) in enumerate(self.pts)]
        if len(dseq) == 1:
            self.isel = 0
        elif dseq[-1] > dseq[-2]:  # guess whether it's an inner or outer probe
            self.isel = min(dseq)[1]
        else:
            self.isel = max(dseq)[1]
        self.isel = 0  # SET BACK TO FIRST POINT OF CONTACT
        if len(self.pts) > 1:
            self.lng = sum(sqrt((self.pts[i][0]-self.pts[i-1][0])**2 + (self.pts[i][1]-self.pts[i-1][1])**2)  for i in range(1, len(self.pts)))
        else:
            self.lng = 0
        self.vel = self.lng/len(self.pts)
    def getselj0(self):
        return float(self.samplelines[self.isel][0])
    def getselj1(self):
        return float(self.samplelines[self.isel][1])
        
        
print("opening file '%s'" % fname)
samplelines = [ (i, ln.split())  for i, ln in enumerate(open(fname))  if i != 0 ]
contactsequences = [ None ]
gapsequence = [ ]
for i, sampleline in samplelines:
    if len(sampleline) == 4:
        if sampleline[3] == '0':
            if contactsequences[-1] is None and len(gapsequence) < Nmergablegap:
                contactsequences.pop()
                contactsequences[-1].ninternalgap += len(gapsequence)
                contactsequences[-1].samplelines.extend(gapsequence)
            if contactsequences[-1] is None:
                contactsequences[-1] = ContactSequence(i, sampleline, len(gapsequence), len(contactsequences)-1)
            else:
                contactsequences[-1].samplelines.append(sampleline)
            gapsequence = [ ]
        else:
            if contactsequences[-1] is not None:
                contactsequences.append(None)
            gapsequence.append(sampleline)


if contactsequences[-1] is None:
    contactsequences.pop()

timegaps = [cs1.i - cs0.i  for cs0, cs1 in zip(contactsequences, contactsequences[1:])]
avgtimegap = sum(timegaps)/len(timegaps)
print("Found %d lines equating to %f minutes" % (len(samplelines), len(samplelines)/halsampleHz/60))
print("Contacts: %d average %f seconds apart, max internalgap %d" % \
      (len(contactsequences), avgtimegap/halsampleHz, 
       max(contactsequence.ninternalgap  for contactsequence in contactsequences)))
del samplelines

# allocate the contact sequences to their disks
[ contactsequence.ApplyConversion(X0)  for contactsequence in contactsequences ]
midx = sum(cs.pts[0][0]  for cs in contactsequences)/len(contactsequences)
midy = sum(cs.pts[0][1]  for cs in contactsequences)/len(contactsequences)
[ contactsequence.guessisel(midx, midy)  for contactsequence in contactsequences ]
#print("iseq of contact sequence", [ contactsequence.isel  for contactsequence in contactsequences ])
#for cs in contactsequences:
#    print(cs.isel, len(cs.pts), [sqrt((x - midx)**2 + (y - midy)**2)  for x, y in cs.pts])

print("Thin out the contact sequences of the double hits")    
print("should project back to find the length of the trajectory before contact to verify")
print("for sure we have a short resample")
print([cs.N  for cs in contactsequences])
print([cs.lng  for cs in contactsequences])
contactsequences = contactsequences[1::2]


for cs in contactsequences:  
    cs.dmid = sqrt((cs.pts[cs.isel][0] - midx)**2 + (cs.pts[cs.isel][1] - midy)**2)
dmids = [ cs.dmid  for cs in contactsequences ]
dmids.sort()
dgaps = [ (dm0 + dm1)/2  for dm0, dm1 in zip(dmids, dmids[1:])  if dm1 - dm0 > 10 ]
print(dgaps)
contactsequencedisks = [[] for i in range(len(dgaps)+1)]
for cs in contactsequences:
    i = 0
    while i < len(dgaps) and cs.dmid > dgaps[i]:
        i += 1
    contactsequencedisks[i].append(cs)

def radsd(C, pts):
    cx, cy = C
    n = len(pts)
    dsqs = [ (cx - x)**2 + (cy - y)**2  for x, y in pts ]
    sumdls = sum(sqrt(dsq)  for dsq in dsqs)
    sumdsqs = sum(dsqs)
    return sumdls/n, sqrt(max(0, sumdsqs*n - sumdls**2))/(n-1)

def radrg(C, pts):
    cx, cy = C
    dls = [ sqrt((cx - x)**2 + (cy - y)**2)  for x, y in pts ]
    return min(dls), max(dls)

def BestCentre(pts):
    cx0 = sum(x  for x, y in pts)/len(pts)
    cy0 = sum(y  for x, y in pts)/len(pts)
    def fun(C):
        return radsd(C, pts)[1]
    b = 5
    bnds = ((cx0-b,cx0+b), (cy0-b,cy0+b))
    res = scipy.optimize.minimize(fun=fun, x0=(cx0, cy0), bounds=bnds)
    return tuple(res.x)

def BestCentreR2(pts, R2):
    cx0 = sum(x  for x, y in pts)/len(pts)
    cy0 = sum(y  for x, y in pts)/len(pts)
    def fun(C):
        cx, cy = C
        dss = [ sqrt((cx - x)**2 + (cy - y)**2)  for x, y in pts ]
        return sum((d - R2)**2  for d in dss)/len(pts)
    b = 5
    bnds = ((cx0-b,cx0+b), (cy0-b,cy0+b))
    res = scipy.optimize.minimize(fun=fun, x0=(cx0, cy0), bounds=bnds)
    return tuple(res.x)


def ApplyRadialFactor(cx, cy, r, p, radialfactor):
    vx, vy = p[0] - cx, p[1] - cy
    vsq = vx**2 + vy**2
    vlen = sqrt(vsq)
    vfac = (vlen - r)*radialfactor / r
    return p[0] + vx*vfac, p[1] + vy*vfac
        
def PlotCircles(X, radialfactor, diskcen=(0,0), matnumber=0):
    dx, dy = diskcen
    for i in range(len(contactsequencedisks)):
        pts = [ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in contactsequencedisks[i]]
        print("Disk %d points %d" % (i, len(contactsequencedisks[i])), end=" ")
        sendactivity("points", points=[(x-dx, y-dy)  for x, y in pts], materialnumber=(matnumber+i)%4)
        cx, cy = BestCentre(pts)
        sendactivity("points", points=[(cx-dx, cy-dy)], materialnumber=(matnumber+i)%4)
        rC, sdC = radsd((cx, cy), pts)
        rmin, rmax = radrg((cx, cy), pts)
        print("Centre (%f,%f) minrad %f maxrad %f" % (cx, cy, rmin, rmax))
        sendactivity("contours", contours=[[(cx+rC*sin(radians(d/5))-dx, cy+rC*cos(radians(d/5))-dy) for d in range(0,360*5+1)]], materialnumber=matnumber)
        #if radialfactor:
        #    sendactivity("contours", contours=[ [ ApplyRadialFactor(cx, cy, rC, procforward(float(sampleline[0]), float(sampleline[1]), X), radialfactor)  for sampleline in cs.samplelines ]  for cs in contactsequencedisks[i] ])
            


sendactivity("clearalltriangles")
sendactivity("clearallpoints")
sendactivity("clearallcontours")
sendactivity("points", points=[cs.pts[0]  for cs in contactsequences])
PlotCircles(X0, 100)


# now minimize on the main set of parameters
def fun(X):
    ptslist = [ [ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in csd ]  for csd in contactsequencedisks ]
    diskcens = [ BestCentre(pts)  for pts in ptslist ]
    rcsdlist = [ radsd(cen, pts)  for cen, pts in zip(diskcens, ptslist) ]
    
    return sum(sd**2  for r, sd in rcsdlist)


# now minimize on the main set of parameters
R2 = 132.525
def fun(X):
    ptsouter = [ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in contactsequencedisks[2] ]
    diskcen = BestCentreR2(ptsouter, R2)
    cx, cy = diskcen
    n = len(ptsouter)
    distss = [ sqrt((cx - x)**2 + (cy - y)**2)  for x, y in ptsouter ]
    return sum((d - R2)**2  for d in distss)/n

print("Initial fun to minimize value", fun(X0))

#res = scipy.optimize.minimize(fun=fun, x0=x0, method='Powell', options={"xtol":0.00000001, "ftol":0.00000001})
tstart = time.time()
bnds = [ (x-5,x+5)  for x in X0 ]



method, matnumber = 'Nelder-Mead', 3
method, matnumber = 'Powell', 1

# 'Nelder-Mead' - 'Powell' - 'CG' - 'BFGS' - 'L-BFGS-B' - 'TNC' - 'COBYLA' - 'SLSQP' - 'dogleg' - 'trust-ncg'
res = scipy.optimize.minimize(fun=fun, x0=X0, bounds=bnds, method=method, options={"xtol":0.00000001, "ftol":0.00000001})
print("Elapsed seconds", time.time() - tstart)
print(res)
Xopt = tuple(res.x)
X = Xopt
diskcen = BestCentreR2([ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in contactsequencedisks[2] ], R2)

print("error %f for %s" % (fun(X0), repr(X0)))
print("error %f for %s" % (fun(Xopt), repr(Xopt)))
PlotCircles(Xopt, 1, diskcen, matnumber)


print("""
void setupconstants(struct arcdata_data *hd) 
{
    hd->a = %f; 
    hd->b = %f; 
    hd->c = %f; 
    hd->d = %f; 
    hd->e = %f; 
    hd->f = %f; 
    hd->ah = %f; 
    hd->tblx = %f; 
    hd->armx = %f; 
}""" % Xopt)


j0s = [cs.getselj0()  for cs in contactsequencedisks[2]]
j1s = [cs.getselj1()  for cs in contactsequencedisks[2]]
j0min, j0max = min(j0s), max(j0s)
j1min, j1max = min(j1s), max(j1s)
nstripes = 10
X = Xopt
sendactivity("contours", contours=[[ procforward(j0min + (j/nstripes)*(j0max - j0min), j1min + (i/100)*(j1max - j1min), X, diskcen)  for i in range(100) ]  for j in range(nstripes)], materialnumber=matnumber)
sendactivity("contours", contours=[[ procforward(j0min + (j/100)*(j0max - j0min), j1min + (i/nstripes)*(j1max - j1min), X, diskcen)  for j in range(100) ]  for i in range(nstripes)], materialnumber=matnumber)

