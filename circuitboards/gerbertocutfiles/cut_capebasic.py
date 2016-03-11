import sys, os  # has to be Python27
sys.path.append("/home/goatchurch/geom3d/circuitcuttingtools")
sys.path.append("/home/goatchurch/geom3d/flatcam")
sys.path.append("/home/goatchurch/geom3d/barmesh")
from gerbergetbits import GetIsolationCuts, GetHoleCuts, GetEdgeCuts
from nongerberbits import ProbingZ, post, sendgcodetotwisted, ReorderConts
from basicgeo import I1

dr = "/home/goatchurch/geom3d/maghullmachinecalibrations/circuitboards/bbcape3"
f = os.path.join(dr, "capebasic_without_gnds-B_Cu.gbl")
drl = os.path.join(dr, "capebasic.drl")
ed = os.path.join(dr, "capebasic_without_gnds-Edge_Cuts.gbr")
prb = os.path.join(dr, "probing6.txt")
ngcF = os.path.join(dr, "bbcape.ngc")
twistedipnumber = "192.168.0.45"


econt = GetEdgeCuts(ed, 1, 0, 0)[0][0]
xlo, xhi = min(p[0]  for p in econt), max(p[0]  for p in econt)
ylo, yhi = min(p[1]  for p in econt), max(p[1]  for p in econt)
print("xyrange", (xlo, xhi), (xhi, yhi))
print("width height", xhi-xlo, yhi-ylo)
dx, dy = -xlo, -ylo
ymid = (yhi+ylo)/2+dy

# scp machinekit@192.168.0.45:machinekit/probing6.txt .
prbz = ProbingZ(prb)
sendactivity(points=[(x,y,z*50)  for (x,y), z in prbz.probeptsmap.items()])
xrg, yrg = I1(prbz.parx.vs[0], prbz.parx.vs[-1]), I1(prbz.pary.vs[0], prbz.pary.vs[-1])
pp = [(xrg.Along(i/50.0), yrg.Along(j/50.0))  for i in range(51)  for j in range(51)]
sendactivity(points=[(x,y,prbz.InterpZ(x,y)*50)  for x,y in pp], materialnumber=1)


fconts, fcontoffsets = GetIsolationCuts(f, 0.22, dx, dy)
sendactivity(contours=fconts, materialnumber=0)
sendactivity(contours=fcontoffsets, materialnumber=1)

holediampts = GetHoleCuts(drl, dx, dy)
sendactivity(points=sum(holediampts.values(), []))

econts, econtoffsets = GetEdgeCuts(ed, 0.4, dx, dy)
econt, econtoffset = econts[0], econtoffsets[0]
sendactivity(contours=[econt], materialnumber=2)
sendactivity(contours=[econtoffset], materialnumber=3)


# copper face isolation tracks
lconts = [[(p[0],ymid*2-p[1])  for p in cont]  for cont in fcontoffsets]
cconts = ReorderConts(lconts)   # this does the 1.8mm runoff at end
post(cconts, -0.1, 2, ngcF, prbz)
sendgcodetotwisted(ngcF, twistedipnumber)

# copper face drill pecks
dconts = [[(p[0],ymid*2-p[1])]  for p in sum(holediampts.values(), [])]
cconts = ReorderConts(dconts)
post(cconts, -0.166, 2, ngcF, prbz)

# copper face edge scoring
eecontoffset = [(p[0],ymid*2-p[1])  for p in econtoffset]
post([eecontoffset], -0.6, 2, ngcF, prbz)


# copper face full depth drills
post(cconts, -3.1, 2, ngcF, prbz)




2.342162
x1,y1 = (89.349, -66.391)
x0,y0 = (0, -40.146)
vx,vy = x1-x0, y1-y0
import math
k = math.atan2(vy, vx)

2.342162-abs(k)+math.pi/2   -> the rotval










