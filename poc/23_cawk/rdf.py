#/bin/sh

python /dev/stdin "$@" <<'EOF'
import sys
import numpy as np
def Read(p):
    return np.genfromtxt(p, delimiter=',', names=True)
def msg(s):
  sys.stderr.write(str(s) + "\n")

rmin = 1.0/192
rmax = 0.1
bins = 50
sval = np.zeros(bins)

for f in sys.argv[1:]:
    msg(f)
    d = Read(f)
    d = d[np.where( (d["r"] < rmax) & (d["r"] > rmin) )]
    x = d["x"]
    y = d["y"]
    z = d["z"]
    r = d["r"]
    n = x.size
    rng = (rmin, rmax)
    for i in range(n):
        dx = x - x[i]
        dy = y - y[i]
        dz = z - z[i]
        dr = np.array([dx, dy, dz])**2
        dr = np.mean(dr, axis = 0)**0.5
        dr = dr[np.where(dr > rmin)]
        val, edge = np.histogram(dr, bins=bins, weights = 1/dr**2, range=rng)
        sval += val

sval /= sval.mean()        
left = edge[:-1]
right = edge[1:]
center = (left + right) * 0.5
np.savetxt(sys.stdout, np.transpose([center, sval, left, right]))
    
EOF
