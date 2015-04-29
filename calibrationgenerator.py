import os

try:  sendactivity
except NameError:
    def sendactivity(*args, **kwargs):  pass
    
sendactivity("clearalltriangles")
sendactivity("clearallpoints")
sendactivity("clearallcontours")

halsampleHz = 1000
fname = "/home/goatchurch/geom3d/maghullmachinecalibrations/data/2015-04-fastprobing.txt"
print("opening file '%s'" % fname)
samplelines = [ (i, ln.split())  for i, ln in enumerate(open(fname))  if i != 0 ]
contactsequences = [ [ ] ]
for i, ls in samplelines:  # maybe file should include line number to make this simpler
    if len(ls) == 4:
        if ls[3] == '0':
            contactsequences[-1].append((i, ls))
        elif contactsequences[-1]:
            contactsequences.append([])
if not contactsequences[-1]:
    contactsequences.pop()
timegaps = [cs1[0][0] - cs0[0][0]  for cs0, cs1 in zip(contactsequences, contactsequences[1:])]
avgtimegap = sum(timegaps)/len(timegaps)

print("Found %d lines equating to %f minutes with %d samples on average %f seconds apart" % \
      (len(samplelines), len(samplelines)/halsampleHz/60, len(contactsequences), avgtimegap/halsampleHz))
           

plotscale = 0.1
cont = [ (float(ls[0])*0.1, float(ls[1])*0.1)  for ls in slns[::10] ]
sendactivity("contours", contours=[cont], materialnumber=1)

pts = [ (float(ls[0])*0.1, float(ls[1])*0.1)  for ls in slns  if ls[3] == '0' ]
sendactivity("points", points=pts, materialnumber=1)

tabarm = [ (float(ls[0]), float(ls[1]))  for ls in slns[::10] ]
tabarm[100]


#pts = [ (float(ls[0])*0.1, float(ls[1])*0.1)  for ls in slns  if ls[3] == '0' ]
#len(pts)
#sendactivity("points", points=pts, materialnumber=1)


a, b, c, d, e, f, h = 745.0, 746.5, 743.5, 695.0, 693.0, 582.625
armx, tblx = 121.5, 221.01165853896501, 248.74565762513299
tablerotoffsetrad = 0.4837

from math import sqrt as SQR
from math import atan as ATN
from math import cos as COS
from math import sin as SIN
from math import acos as ACOS

def ProcForward(a, b, c, d, e, f, h, tbl, arm): # Screwlength to Vertex coords
    #DIM AS DOUBLE g,aH,aX,aF,aG,aTBL,aARM'variables
    aB=ACOS(((a**2)+(c**2)-(b**2))/(2*a*c))
    aC=ACOS(((a**2)+(b**2)-(c**2))/(2*a*b))
    aH=ACOS(((f*f)+(e*e)-(h*h))/(2*f*e))
    aTBL=ACOS(((d*d)+(b*b)-(tbl*tbl))/(2*d*b))     #<<<<<<< tbl
    aARM=ACOS(((c*c)+(e*e)-(arm*arm))/(2*c*e))     #<<<<<<< arm
    aG  =aB-aARM+aH
    aF  =ATN((f*(SIN(aG)))/(a-(f*(COS(aG)))))
    aX  =aTBL-(aC-aF)-tablerotoffsetrad
    g   =SQR((f*f)+(a*a)-(2*f*a*(COS (aG))))
    x   =(g*SIN(aX))+110                           #x >>>>>>>
    y   =572-(g*COS(aX))                           #y >>>>>>>
    
    return x, y

css = [ ProcForward(a, b, c, d, e, f, h, tbl, arm)  for tbl, arm in tabarm ]
len(css)
css[1000]
sendactivity("contours", contours=[[(x*0.1, y*0.1)  for x, y in css]], materialnumber=3)



