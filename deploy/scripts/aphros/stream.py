try:
    import numpy as np
except ImportError:
    pass

try:
    import scipy
    import scipy.interpolate
    import scipy.sparse as sp
    import scipy.sparse.linalg
except ImportError:
    pass

verbose = False


def Log(msg):
    if verbose:
        print(str(msg))


def report(re, n="residual"):
    re = re.flatten()
    nr = np.linalg.norm
    l = len(re)
    Log("{:}: L1: {:}, L2: {:}, Linf: {:}".format(n,
                                                  nr(re, ord=1) / l,
                                                  nr(re, ord=2) / (l**0.5),
                                                  nr(re, ord=np.inf)))


# vorticity: dv/dx - du/dy
# ADHOC: periodic in x
def vort(u, v):
    vx = 0.5 * (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1))
    uy = 0.5 * (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0))
    om = vx - uy
    #om[:,0] = 0.
    #om[:,-1] = 0.
    om[0, :] = 0.
    om[-1, :] = 0.
    return om


# solve Posson equation for
# laplace p = -vort(u,v)
# stream function p (psi) defined via u = dp/dy, v = -dp/dx
def stream_jacobi(u, v, itmax=50000, tol=1e-8):
    Log("stream_jacobi")
    # indexing: u[y,x]
    f = -vort(u, v)
    s = np.zeros_like(f)
    it = 0
    r = tol + 1.
    lastit = 0
    while r > tol and it < itmax:
        pold = s
        sor = 1.0
        resid = (np.roll(s, -1, axis=0) + np.roll(s, 1, axis=0) + np.roll(
            s, -1, axis=1) + np.roll(s, 1, axis=1) - 4. * s - f) / 4.
        s = pold + sor * resid
        # boundary conditions: derivatives
        s[:, 0] = s[:, 1] + (v[:, 0] + v[:, 1]) * 0.5
        s[:, -1] = s[:, -2] - (v[:, -1] + v[:, -2]) * 0.5
        s[0, :] = s[1, :] - (u[0, :] + u[1, :]) * 0.5
        s[-1, :] = s[-2, :] + (u[-1, :] + u[-2, :]) * 0.5
        # center at zero
        s -= s.mean()
        r = (resid[1:-1, 1:-1]**2).mean()**0.5
        # report
        if it >= lastit + itmax // 10:
            Log("i={:06d}, r={:.5e}".format(it, r))
            lastit = it
        it += 1
    return s


# solve Posson equation for
# laplace p = -vort(u,v)
# stream function p (psi) defined via u = dp/dy, v = -dp/dx
# ADHOC: periodic in x
def stream_direct(u, v):
    Log("stream_direct")
    f = -vort(u, v)

    i = np.arange(len(f.flatten())).reshape(f.shape)
    ones = np.ones_like(i).astype(float)

    # inner: laplacian
    dc = -4. * ones.copy()
    dxp = ones.copy()
    dxm = ones.copy()
    dyp = ones.copy()
    dym = ones.copy()

    ic = i
    ixp = np.roll(i, -1, axis=1)
    ixm = np.roll(i, 1, axis=1)
    iyp = np.roll(i, -1, axis=0)
    iym = np.roll(i, 1, axis=0)

    # boundary conditions: derivatives
    # laplacian = (xp - c) - (c - xm) + (yp - c) - (c - ym)
    # subtract gradients at the sides
    #dxp[:,-1] -= 1. ; dc[:,-1] += 1.
    #dxm[:, 0] -= 1. ; dc[:, 0] += 1.
    dyp[-1, :] -= 1.
    dc[-1, :] += 1.
    dym[0, :] -= 1.
    dc[0, :] += 1.

    # add given values
    #f[:,-1] += (v[:,-1] + v[:,-2]) * 0.5
    #f[:, 0] -= (v[:, 0] + v[:, 1]) * 0.5
    f[-1, :] += -(u[-1, :] + u[-2, :]) * 0.5
    f[0, :] -= -(u[0, :] + u[1, :]) * 0.5

    # combine
    ld = (dc, dxp, dxm, dyp, dym)
    data = np.stack(ld, axis=-1)
    li = (ic, ixp, ixm, iyp, iym)
    indices = np.stack(li, axis=-1)
    rowsize = len(li) * np.ones_like(ic)

    # flatten
    data = data.flatten()
    indices = indices.flatten()
    rowsize = rowsize.flatten()
    indptr = np.concatenate(([0], np.cumsum(rowsize)))

    # fix s[0] = 0
    data[0] += 1.

    # solve ax=b
    a = sp.csr_matrix((data, indices, indptr))
    b = f.flatten()
    s = sp.linalg.spsolve(a, b)
    report(a.dot(s) - b)
    s = s.reshape(f.shape)

    # compute u,v
    su = 0.5 * (np.roll(s, -1, axis=0) - np.roll(s, 1, axis=0))
    sv = -0.5 * (np.roll(s, -1, axis=1) - np.roll(s, 1, axis=1))
    su[0, :] = (s[1, :] - s[0, :])
    su[-1, :] = (s[-1, :] - s[-2, :])
    report((su - u) / (abs(u)).mean(), "|u-su|/mean(|u|)")
    report((sv - v) / (abs(v)).mean(), "|v-sv|/mean(|v|)")

    return s


stream = stream_direct
#stream = stream_jacobi


def interp_to_uniform(xp, yp, fields=[], n=None):
    if n is None:
        n = int(len(xp)**0.5 + 0.5)

    x1 = np.linspace(xp.min(), xp.max(), n)
    y1 = np.linspace(yp.min(), yp.max() - 1e-6 * yp.ptp(), n)
    # ADHOC: avoid fill_value in interpolate

    x, y = np.meshgrid(x1, y1)
    pts = np.vstack((xp, yp)).T

    r = []
    for fp in fields:
        f = scipy.interpolate.griddata(pts, fp, (x, y), method="cubic")
        assert np.all(np.isfinite(f))
        r.append(f)

    return x1, y1, r
