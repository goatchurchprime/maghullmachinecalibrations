j0rg = [100, -315]
j1rg = [-218, 100]

def Along(lam, rg):
    return rg[0]*(1-lam) + rg[1]*lam

lams = [ 0.05+0.9*i/3 for i in range(4) ]


pts = [ (Along(l0, j0rg), Along(l1, j1rg)) for l0 in lams  for l1 in lams ]

zsafe = 50
ztopplus = 22
ztop = 20
zdepths = [5, 10, 15]
feed = 200

pt = pts[0]

def PeckCycle(pt):
    res = [ ]
    res.append("G0Z%d" % zsafe)
    res.append("X%0.3fY%0.3f" % pt)
    res.append("G1Z%dF%d" % (ztopplus, feed))
    for zdepth in reversed(zdepths):
        res.append("Z%d" % (ztop))
        res.append("Z%d" % (zdepth))
    res.append("")
    return res

fout = open("/home/goatchurch/geom3d/maghullmachinecalibrations/joint4.ngc", "w")
for pt in pts:
    fout.write("\n".join(PeckCycle(pt)))
fout.write("M2\n")
fout.close()

