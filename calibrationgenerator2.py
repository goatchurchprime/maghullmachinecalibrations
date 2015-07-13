from math import sqrt, atan, cos, sin, acos, radians
import os, time, scipy.optimize

# circle sizes = 132.525, 79.69, 17.055

try:  sendactivity
except NameError:
    def sendactivity(*args, **kwargs):  pass
    
fname = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/2015-05-22-diskprobing.txt"
fnamelin = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/2015-06-01-linearprobe.txt"

X0 = (744.57703517328832, 746.19009986079459, 741.95547170764314, 695.42905812855656, 693.07538973655983, 582.64403220610257, 0.0764480526, 690, 400)
#X0 = (743.3133003594869, 756.82976800937217, 743.43834112653599, 698.16794077213683, 688.61101385919824, 580.7463415098855, 0.097855119579137784, 518.65650220023929, 408.6731700849632)

halsampleHz = 1000
Nmergablegap = 160  # when finding samples in contact
aXoffset = 0#-0.4837

def GetContactSequence(fname):
    print("opening file '%s'" % fname)
    samplelines = [ sl  for sl in (ln.split()  for ln in open(fname))  if len(sl) == 4 ]
    print("Found %d lines equating to %f minutes" % (len(samplelines), len(samplelines)/halsampleHz/60))
    contactsequences = [ ]
    for i, sampleline in enumerate(samplelines):
        if sampleline[3] == '0':
            if contactsequences and i < contactsequences[-1][1] + Nmergablegap:
                contactsequences[-1][1] = i
            else:
                contactsequences.append([i, i])
    print("Contacts: %d" % (len(contactsequences)))      
    return samplelines, contactsequences

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


def AllocatePointsToDisksIndexes(contactpts):
    midx = sum(pt[0]  for pt in contactpts)/len(contactpts)
    midy = sum(pt[1]  for pt in contactpts)/len(contactpts)
    dmids = [ (sqrt((pt[0] - midx)**2 + (pt[1] - midy)**2), i)  for i, pt in enumerate(contactpts) ]
    dmids.sort()
    dgaps = [ (dm0 + dm1)/2  for (dm0, i0), (dm1, i1) in zip(dmids, dmids[1:])  if dm1 - dm0 > 10 ]
    # circle sizes = 132.525, 79.69, 17.055
    diskindexes = [ []  for i in range(len(dgaps)+1) ]
    for dmid, i in dmids:
        j = 0
        while j < len(dgaps) and dmid > dgaps[j]:
            j += 1
        diskindexes[j].append(i)
    return diskindexes
    
def BestCentreR2(pts, R2):
    cx0, cy0 = sum(x  for x, y in pts)/len(pts), sum(y  for x, y in pts)/len(pts)
    def centrequality(c):
        dss = [ sqrt((c[0] - x)**2 + (c[1] - y)**2)  for x, y in pts ]
        return sum((d - R2)**2  for d in dss)/len(pts)
    b = 5
    bnds = ((cx0-b,cx0+b), (cy0-b,cy0+b))
    res = scipy.optimize.minimize(fun=centrequality, x0=(cx0, cy0), bounds=bnds)
    return tuple(res.x), centrequality(res.x)

def BestLineNormal(pts):
    sx = sum(x  for x, y in pts)
    sy = sum(y  for x, y in pts)
    sx2 = sum(x**2  for x, y in pts)
    sxy = sum(x*y  for x, y in pts)
    m = (sxy*len(pts) - sx*sy) / (sx2*len(pts) - sx**2)
    c = (sy - m*sx)/len(pts)
    nveclen = sqrt(m**2 + 1)
    return (m/nveclen, -1/nveclen)

fac745 = 0.1

def OptimizeXagainstR2(X0, outerdisksamples, R2, linesamples, method):
    tstart = time.time()
    bnds = [ (x-5,x+5)  for x in X0 ]
    def fun1(X):
        ptsouter = [ procforward(j0, j1, X)  for j0, j1 in outerdisksamples ]
        return BestCentreR2(ptsouter, R2)
    def fun(X):
        var = 0
        if linearsamples:
            linepts = [ procforward(j0, j1, X)  for j0, j1 in linearsamples ]
            nvec = BestLineNormal(linepts)
            ds = [ nvec[0]*x + nvec[1]*y  for x, y in linepts ]
            var = (sum(d**2  for d in ds) - sum(ds)**2/len(ds))/len(ds)
        var745 = (X[0]-745)**2 + (X[1]-745)**2 + (X[2]-745)**2
        return fun1(X)[1] + var + var745*fac745
    res = scipy.optimize.minimize(fun=fun, x0=X0, bounds=bnds, method=method, options={"xtol":0.00000001, "ftol":0.00000001, "Dmaxiter":9000})
    print("Elapsed seconds", time.time() - tstart)
    print(res)
    Xopt = tuple(res.x)
    return Xopt, fun1(Xopt)

    
# this does the plotting and measuring after optimization stuff
def OptMeth(X0, meth):
    k0 = OptimizeXagainstR2(X0, outerdisksamples, R2, linearsamples, meth)
    X = k0[0]
    ptsouter0 = [ procforward(j0, j1, X)  for j0, j1 in outerdisksamples ]
    c = BestCentreR2(ptsouter0, R2)[0]
    ds = [ sqrt((p[0] - c[0])**2 + (p[1] - c[1])**2)  for p in ptsouter0 ]
    print("Centre (%f,%f) minrad %f maxrad %f" % (c[0], c[1], min(ds), max(ds)))
    matnumber = len(meth)%4
    sendactivity("contours", contours=[[(R2*sin(radians(d/5)), R2*cos(radians(d/5))) for d in range(0,360*5+1)]], materialnumber=matnumber)
    sendactivity("points", points=[(x-c[0], y-c[1])  for x, y in ptsouter0], materialnumber=matnumber)
    
    linepts = [ procforward(j0, j1, X)  for j0, j1 in linearsamples ]
    nvec = BestLineNormal(linepts)
    nvecdots = [ nvec[0]*x + nvec[1]*y  for x, y in linepts ]
    print("linear width is ", max(nvecdots) - min(nvecdots))
    sendactivity("points", points=[(x-c[0], y-c[1])  for x, y in linepts], materialnumber=matnumber)
    nvc = sum(nvecdots)/len(nvecdots)
    sendactivity("contours", contours=[[(nvec[0]*nvc + nvec[1]*d - c[0], nvec[1]*nvc - nvec[0]*d - c[1])  for d in [-1000,1000]]], materialnumber=matnumber)
    
    # plot the constant joint lines
    j0s = [j0  for j0, j1 in outerdisksamples]
    j1s = [j1  for j0, j1 in outerdisksamples]
    j0min, j0max = min(j0s), max(j0s)
    j1min, j1max = min(j1s), max(j1s)
    nstripes = 10
    sendactivity("contours", contours=[[(p[0]-c[0], p[1]-c[1])  for p in [ procforward(j0min + (j/nstripes)*(j0max - j0min), j1min + (i/100)*(j1max - j1min), X)  for i in range(100) ]]  for j in range(nstripes)], materialnumber=matnumber)
    sendactivity("contours", contours=[[(p[0]-c[0], p[1]-c[1])  for p in [ procforward(j0min + (j/100)*(j0max - j0min), j1min + (i/nstripes)*(j1max - j1min), X)  for j in range(100) ]]  for i in range(nstripes)], materialnumber=matnumber)
    
    return X
    
    
# main code

# parse out the contacts from the halsampler file
samplelines, contactsequences = GetContactSequence(fname)
contactsequences = contactsequences[1::2]   # thin out double values

contactpts = [ procforward(float(samplelines[i0][0]), float(samplelines[i0][1]), X0)  for i0, i1 in contactsequences ]

# select the outer disk (and the outerdisksamples containing the joint values)
sendactivity("points", points=contactpts)
diskindexes = AllocatePointsToDisksIndexes(contactpts)
outerdiskpts = [contactpts[i]  for i in diskindexes[2]]
sendactivity("points", points=outerdiskpts, materialnumber=1)
outerdisksamples = [ (float(samplelines[contactsequences[i][0]][0]), float(samplelines[contactsequences[i][0]][1]))  for i in diskindexes[2] ]
sendactivity("points", points=outerdisksamples, materialnumber=2)

# fit the outer disk with current estimate
R2 = 132.525
c, q = BestCentreR2(outerdiskpts, R2)
print("minimizing from", q)
sendactivity("contours", contours=[[(c[0]+R2*sin(radians(d/5)), c[1]+R2*cos(radians(d/5))) for d in range(0,360*5+1)]], materialnumber=1)

# read the straight rule line sampling file, to add into the optimization function
linsamplelines, lincontactsequences = GetContactSequence(fnamelin)
linearsamples = [ (float(linsamplelines[i0][0]), float(linsamplelines[i0][1]))  for i0, i1 in lincontactsequences ]

ffff
print("Introducing fac745 =", fac745, "to force X[0:3] to be 745")
sendactivity("clearalltriangles")
sendactivity("clearallpoints")
sendactivity("clearallcontours")
X = OptMeth(X0, 'Nelder-Mead')
print(X)
X = OptMeth(X0, 'Powell')
print(X)
    
print("void setupconstants(struct arcdata_data *hd)\n{\r\n%s\r\n}" % "\r\n".join("    hd->%s = %f;" % (s, v)  for s, v in zip(("a", "b", "c", "d", "e", "f", "ah", "tblx", "armx"), X)))
