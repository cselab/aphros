import os


def GetVal(path, colname, t0):
    s = os.path.join(path, "stat.dat")

    res = 0.

    # best
    tb = None
    vb = None

    with open(s) as f:
        head = f.readline().split()
        it = head.index('t')
        iv = head.index(colname)
        for l in f:
            t = float(l.split()[it])
            v = l.split()[iv]
            if tb is None or abs(t - t0) < abs(tb - t0):
                tb = t
                vb = v

    return float(vb)
