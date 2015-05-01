from math import sqrt, atan, cos, sin, acos, radians
import scipy.optimize
import os, time

try:  sendactivity
except NameError:
    def sendactivity(*args, **kwargs):  pass
    
fname = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/2015-04-slowprobing.txt"
#fname = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/2015-04-fastprobing.txt"

X0 = (745.0, 746.5, 743.5, 695.0, 693.0, 582.625, 121.493035)  # initial guess
X0 = (743.60582782689789, 746.50166754871304, 743.46653193215002, 693.44977777514146, 693.54373345096121, 582.74575870510228, 121.21098136397796)
X0 = (743.60427577698113, 746.49690047566355, 743.46654290854087, 693.46920163003449, 693.54373390841067, 582.74599745793853, 121.21185734841112)

halsampleHz = 1000
Nmergablegap = 160  # when finding samples in contact

def procforward(j0, j1, X):
    a, b, c, d, e, f, h = X
    aB = acos(((a**2)+(c**2)-(b**2))/(2*a*c))
    aC = acos(((a**2)+(b**2)-(c**2))/(2*a*b))
    aH = acos(((f*f)+(e*e)-(h*h))/(2*f*e))
    aTBL = acos(((d*d)+(b*b)-(j0*j0))/(2*d*b))
    aARM = acos(((c*c)+(e*e)-(j1*j1))/(2*c*e))
    aG = aB-aARM+aH
    aF = atan((f*(sin(aG)))/(a-(f*(cos(aG)))))
    aX = aTBL-(aC-aF)
    g = sqrt((f*f)+(a*a)-(2*f*a*(cos(aG))))
    x = g*sin(aX)
    y = -g*cos(aX)
    return x, y


class ContactSequence:
    def __init__(self, i, sampleline, nprevgap):
        self.i = i
        self.samplelines = [ sampleline ]
        self.nprevgap = nprevgap
        self.ninternalgap = 0
    def ApplyConversion(self, X):
        self.pts = [ procforward(float(sampleline[0]), float(sampleline[1]), X)  for sampleline in self.samplelines ]
        self.isel = int((len(self.samplelines)+0.5)/2)  # selected point to use from the direction
    def guessisel(self, mx, my):
        dseq = [((x - mx)**2 + (y - my)**2, i)  for i, (x, y) in enumerate(self.pts)]
        if dseq[-1] > dseq[-2]:  # guess whether it's an inner or outer probe
            self.isel = min(dseq)[1]
        else:
            self.isel = max(dseq)[1]
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
                contactsequences[-1] = ContactSequence(i, sampleline, len(gapsequence))
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
midx = sum(cs.pts[cs.isel][0]  for cs in contactsequences)/len(contactsequences)
midy = sum(cs.pts[cs.isel][1]  for cs in contactsequences)/len(contactsequences)
print([ contactsequence.isel  for contactsequence in contactsequences ])
[ contactsequence.guessisel(midx, midy)  for contactsequence in contactsequences ]
print([ contactsequence.isel  for contactsequence in contactsequences ])
#for cs in contactsequences:
#    print(cs.isel, len(cs.pts), [sqrt((x - midx)**2 + (y - midy)**2)  for x, y in cs.pts])
    


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
        
def PlotCircles(X):
    for i in range(len(contactsequencedisks)):
        pts = [ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in contactsequencedisks[i]]
        print("Disk %d points %d" % (i, len(contactsequencedisks[i])), end=" ")
        sendactivity("points", points=pts, materialnumber=i)
        cx, cy = BestCentre(pts)
        sendactivity("points", points=[(cx, cy)], materialnumber=i)
        rC, sdC = radsd((cx, cy), pts)
        rmin, rmax = radrg((cx, cy), pts)
        print("Centre (%f,%f) minrad %f maxrad %f" % (cx, cy, rmin, rmax))
        sendactivity("contours", contours=[[(cx+rC*sin(radians(d/5)), cy+rC*cos(radians(d/5))) for d in range(0,360*5+1)]], materialnumber=1)

# plot with current parameters    
sendactivity("clearalltriangles")
sendactivity("clearallpoints")
sendactivity("clearallcontours")
sendactivity("contours", contours=[cs.pts  for cs in contactsequences])
PlotCircles(X0)

# now minimize on the main set of parameters
def fun(X):
    ptslist = [ [ procforward(cs.getselj0(), cs.getselj1(), X)  for cs in csd ]  for csd in contactsequencedisks ]
    diskcens = [ BestCentre(pts)  for pts in ptslist ]
    rcsdlist = [ radsd(cen, pts)  for cen, pts in zip(diskcens, ptslist) ]
    return sum(sd**2  for r, sd in rcsdlist)

print("Initial fun to minimize value", fun(X0))

#res = scipy.optimize.minimize(fun=fun, x0=x0, method='Powell', options={"xtol":0.00000001, "ftol":0.00000001})
tstart = time.time()
bnds = [ (x-5,x+5)  for x in X0 ]
res = scipy.optimize.minimize(fun=fun, x0=X0, bounds=bnds, method='Powell', options={"xtol":0.00000001, "ftol":0.00000001})
print("Elapsed seconds", time.time() - tstart)
print(res)
Xopt = tuple(res.x)

print("error %f for %s" % (fun(X0), repr(X0)))
print("error %f for %s" % (fun(Xopt), repr(Xopt)))
PlotCircles(Xopt)


# numbers from the original mk6skins.c file, with h back-calculated to match aH
#X00 = (744.942619, 746.50577, 743.340584, 695.063413, 692.996544, 582.632544, 121.49303499999954)
#a, b, c, d, e, f, h = X00
a, b, c, d, e, f, h = Xopt

print("Code to paste into mk6skins.c:\n")
print("#define DEFAULT_A %.10f" % a)
print("#define DEFAULT_B %.10f" % b)
print("#define DEFAULT_C %.10f" % c)
print("#define DEFAULT_D %.10f" % d)
print("#define DEFAULT_E %.10f" % e)
print("#define DEFAULT_F %.10f" % f)
aB = acos(((a**2)+(c**2)-(b**2))/(2*a*c))
aC = acos(((a**2)+(b**2)-(c**2))/(2*a*b))
aH = acos(((f*f)+(e*e)-(h*h))/(2*f*e))
print("#define DEFAULT_AB %.10f" % aB)
print("#define DEFAULT_AC %.10f" % aC)
print("#define DEFAULT_AH %.10f" % aH)


