
x, y = 200, -190
x, y = 185, -210
x, y = 120, -40
x, y = 0, 0

zlo, zhi = -8, 1
r = (3.125-3)*0.5+0.015
z = zhi
print("G1Z%.3fF800" % zhi)
print("G1X%.3fY%.3fF800" % (x+r, y))
while z > zlo-1:
    print("G2X%.3fY%.3fZ%.3fI%.3fJ%.3fF200" % (x-r, y, max(z,zlo), -r, 0))
    z -= 0.1
    print("G2X%.3fY%.3fZ%.3fI%.3fJ%.3fF200" % (x+r, y, max(z,zlo), r, 0))
    z -= 0.1
print("G1X%.3fY%.3fF800" % (x, y))
print("G1Z%.3fF800" % zhi)

print("M2")



