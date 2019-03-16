function vof_trapz(M, N, f, param, x, y, u, v) {
    var i, j;
    var p, q, dx, dy;
    var s, a, b;
    if (v === undefined)
        throw new Error("v is undefined")
    if (typeof f !== "function")
        throw new Error("f should be a function:" + f)
    s = 0
    dx = (u - x)/M
    dy = (v - y)/N
    for (i = 0; i < M + 1; i++) {
        a = (i == 0 || i == M) ? 1 : 2
        p = x + i*dx
        for (j = 0; j < N + 1; j++) {
            b = (j == 0 || j == N) ? 1 : 2
            q = y + j*dy
            s += a * b * f(p, q, param)
        }
    }
    return s/(4*N*M)
}


function Vof(dx, dy, f, param) {
    if (typeof f !== "function")
        throw new Error("f should be a function: " + f)
    this.dx = dx
    this.dy = dy
    this.f  = f
    this.param = param
    this.cell = function(i, j) {
        var M = 20, N = 20
        var x, y, u, v
        x = this.dx*i; y = this.dy*j
        u = this.dx*(i + 1); v = this.dy*(j + 1)
        return vof_trapz(M, N, this.f, this.param, x, y, u, v)
    }

    this.grid = function(M, N, /**/ a) {
        var i, j
        if (!Array.isArray(a))
            throw new Error("a is not an array: " + a)
        for (i = 0; i < M; i++)
            for (j = 0; j < N; j++)
                a[i][j] = this.cell(i, j)
    }
}
