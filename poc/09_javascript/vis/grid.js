function Clip(a,l,h) {
    if (isFinite(a)) {
        return Math.min(Math.max(a, l), h);
    }
    return 0.0
}

function sqr(a)  { return a*a }
function abs(a)  { return a > 0 ? a : -a }
function sign(a) { return  a == 0 ? 0 : (a >  0 ? 1 : -1) }

function ScreenToIdx(x, y) {
    var i, j
    i = Clip(Math.floor((x - base.x) / w), 0, nx - 1)
    j = ny - Clip(Math.floor((y - base.y) / w), 0, ny - 1) - 1
    return [i, j]
}

function IsGrid(x, y) {
    var i, j;
    i = Math.floor((x - base.x) / w);
    j = ny - Math.floor((y - base.y) / w) - 1;
    return i >= 0 && i < nx && j >= 0 && j < ny;
}

function IdxToScreen(i, j) {
    var x, y
    x = base.x + i * w
    y = base.y + (ny - j - 1) * w
    return [x, y]
}

function PhysToScreen(xp, yp) {
    var x, y
    x = base.x + xp * w
    y = base.y + (ny - yp) * w
    return [x, y]
}

function ScreenToPhys(x, y) {
    var xp, yp
    xp = (x - base.x) / w
    yp = ny - (y - base.y) / w
    return [xp, yp]
}

function RGBToHex(r,g,b) {
    r = Clip(r, 0, 255);
    g = Clip(g, 0, 255);
    b = Clip(b, 0, 255);

    r = r.toString(16);
    g = g.toString(16);
    b = b.toString(16);

    if (r.length == 1)
        r = "0" + r;
    if (g.length == 1)
        g = "0" + g;
    if (b.length == 1)
        b = "0" + b;
    return "#" + r + g + b;
}

function F(a) {
    //return a * a;
    return a;
}

function Finv(a) {
    //return Math.sqrt(a);
    return a;
}

function IncCell() {
    var wk = wx * kbar

    var ij = ScreenToIdx(start.x, start.y)
    var i = ij[0]
    var j = ij[1]


    var dirx = (abs(mouse.x - start.x) > abs(mouse.y - start.y));
    // displacement [px]
    var dy = (dirx ? mouse.x - start.x : mouse.y - start.y)
    // range of u correction: [-dum,dup]
    var dup = 1.0 - ustart[i][j]
    var dum = ustart[i][j]
    // normalize
    var dyn = dy / wk

    var dyp = Finv(dup) * wk
    var dym = Finv(dum) * wk

    dy = (dirx ? Clip(dy, -dym, dyp) : Clip(dy, -dyp, dym))

    // correction
    var du = F(abs(dyn)) * sign(dyn);

    u[i][j] = Clip(ustart[i][j] + du, 0.0, 1.0);
    return [dy, dym, dyp, dirx]
}


function DrawGridFill(q, x, y) {
    ctx.strokeRect(x, y, w, w);
    q = Math.floor((1-q) * 255 + q * 127);
    ctx.fillStyle = RGBToHex(q, q, q);
    ctx.fillRect(x, y, w, w);
}
function DrawGridText(q, x, y) {
    ctx.font = "20px Arial";
    ctx.fillStyle = "#000000";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(q.toFixed(2), x + w * 0.5, y + w * 0.5);
}
function DrawGrid(u) {
    var i, j, q, xy
    for (j = 0 ; j < ny; j++) {
        for (i = 0 ; i < nx; i++) {
            q = u[i][j];
            xy = IdxToScreen(i, j)
            DrawGridFill(q, xy[0], xy[1])
            /* DrawGridText(q, xy[0], xy[1]) */
        }
    }
}

function DrawLines(ee, u) {
    var i, j;
    for (j = 0 ; j < ny; j++) {
        for (i = 0 ; i < nx; i++) {
            var e = ee[i][j]
            if (e === undefined) continue
            var xya = PhysToScreen(e[0], e[1])
            var xyb = PhysToScreen(e[2], e[3])
            var xa = xya[0], ya = xya[1]
            var xb = xyb[0], yb = xyb[1]

            // unit normal
            var mx = -(yb - ya)
            var my = (xb - xa)
            var mm = Math.sqrt(mx * mx + my * my)
            mx /= mm
            my /= mm
            // shift [px]
            var sh = ctx.lineWidth * 0.2
            var dx = sh * mx
            var dy = sh * my

            ctx.strokeStyle = "#ffffff"
            ctx.beginPath();
            ctx.moveTo(xa - dx, ya - dy)
            ctx.lineTo(xb - dx, yb - dy)
            ctx.closePath()
            ctx.stroke()

            ctx.strokeStyle = "#000000"
            ctx.beginPath();
            ctx.moveTo(xa + dx, ya + dy)
            ctx.lineTo(xb + dx, yb + dy)
            ctx.closePath()
            ctx.stroke()
        }
    }
}

function DrawString(pp, color) {
    ctx.beginPath();
    var i
    for (i = 1 ; i < pp.length; ++i) {
        var xya = PhysToScreen(pp[i - 1][0], pp[i - 1][1])
        var xyb = PhysToScreen(pp[i][0], pp[i][1])
        var xa = xya[0], ya = xya[1]
        var xb = xyb[0], yb = xyb[1]

        ctx.strokeStyle = color;
        ctx.beginPath();
        ctx.moveTo(xa, ya)
        ctx.lineTo(xb, yb)
        ctx.closePath();
        ctx.stroke()
    }

    for (i = 0 ; i < pp.length; ++i) {
        var xy = PhysToScreen(pp[i][0], pp[i][1])
        var x = xy[0], y = xy[1]

        ctx.strokeStyle = "black";
        ctx.beginPath();
        ctx.arc(x, y, 3, 0, 2 * Math.PI);
        ctx.closePath();
        ctx.stroke();
    }
    ctx.closePath();
}

function DrawBar(dy, dym, dyp, dirx) {
    var q = wbar
    var qh = q / 2
    var qq = q / 2
    var qqh = qq / 2
    ctx.fillStyle = '#a0a0a0';
    if (dirx) {
        ctx.fillRect(start.x+dyp, start.y-qh, -dyp-dym-qq, q);
    } else {
        ctx.fillRect(start.x-qh, start.y +dym, q, -dyp-dym-qq);
    }
    ctx.fillStyle = 'blue';
    if (dirx) {
        ctx.fillRect(start.x-qqh, start.y-qqh, dy, qq);
    } else {
        ctx.fillRect(start.x-qqh, start.y-qqh, qq, dy);
    }
    ctx.fillStyle = 'green';
    ctx.fillRect(start.x-qh, start.y-qh, q, q);
    //ctx.fillStyle = 'red';
    //ctx.fillRect(start.x-qh, start.y + dy-qh, q, q);
}

function GridToText(u, ny) {
    var t = ""
    var i, j
    for (j = 0; j < ny; ++j) {
        t += (j == 0 ? "" : "\n")
        for (i = 0; i < nx; ++i) {
            t += (i == 0 ? "" : " ") + u[i][ny - j - 1].toFixed(3)
        }
    }
    return t
}

function TextToGrid(text, nx, ny, /**/ u) {
    var t = text.replace(/\s+/g, " ")
    var v = t.split(" ")
    var m = 0
    var i, j
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            u[i][ny - j - 1] = Clip(parseFloat(v[m++]), 0.0, 1.0)
        }
    }
    return u
}

function UpdateText(t) {
    textarea.value = t
}

function UpdateGrid() {
    var t = textarea.value
    TextToGrid(t, nx, ny, u)
}

function DrawSelect(i, j, current) {
    var xy = IdxToScreen(i, j)
    var x = xy[0] + w/2, y = xy[1] + w/2;
    var q = w / 10
    ctx.fillStyle = (current ? "green" : "red")
    ctx.fillRect(x - q, y - q, 2 * q, 2 * q);
}

function DrawSelectOld(i, j) {
    DrawSelect(i, j, false)
}

function DrawSelectCurrent(i, j) {
    DrawSelect(i, j, true)
}

function DrawSelectPos(xp, yp) {
    var xy = PhysToScreen(xp, yp)
    var x = xy[0], y = xy[1];
    ctx.strokeStyle = "black";
    ctx.beginPath();
    ctx.arc(x, y, circrad * w, 0, 2 * Math.PI);
    ctx.closePath();
    ctx.stroke();
}

function DrawStringRun(i0, j0, color) {
    var AX = 0, AY = 1, BX = 2, BY = 3
    if (ends[i0][j0] !== undefined) {
        var end = partstr_cell_ends(nx, ny, i0, j0, ends)
        var ne = end.length
        var nh = 4
        var hp = 4.0 / (2.0*nh)
        var eta = 0.5
        var eps = 1e-5
        var itermax = 20
        var t = 0.0
        var n = 2*nh + 1
        var partstr = new Partstr(nh, hp, eta)
        var e = end[0]
        var p = [(e[AX] + e[BX]) * 0.5, (e[AY] + e[BY]) * 0.5]
        var dx = e[BX] - e[AX]
        var dy = e[BY] - e[AY]
        var a = Math.atan2(dy, dx)
        partstr.start(ne, end, a, t, p)
        partstr.converge(eps, itermax)
        DrawString(partstr.xx, color)
    }
}


function DrawAll() {
    var AX = 0, AY = 1, BX = 2, BY = 3

    matrix_halo_zero(nx, ny, 2, u)
    ends = matrix_new(nx, ny)
    partstr_vof_ends(nx, ny, u, ends)
    DrawGrid(u)
    DrawSelectOld(i0, j0)
    DrawLines(ends, u)
    DrawStringRun(i0, j0, "red")
}

function Clear() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
}

var onPaint = function() {
    Clear()

    DrawAll()

    if (selectcell) {
        if (IsGrid(mouse.x, mouse.y)) {
            var ij = ScreenToIdx(mouse.x, mouse.y)
            DrawSelectCurrent(ij[0], ij[1])
            DrawStringRun(ij[0], ij[1], "green")
        }
    } else if (selectpos) {
        var xyp = ScreenToPhys(mouse.x, mouse.y)
        DrawSelectPos(xyp[0], xyp[1])
    } else if (mousepressed) { // dragging cell
        var dd = IncCell()
        DrawBar(dd[0], dd[1], dd[2], dd[3])
        UpdateText(GridToText(u, ny))
    }
};

function StartSelect() {
    StopSelect()
    StopSelectPos()

    selectcell = true;
    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
}

function ApplySelect() {
    var x = mouse.x, y = mouse.y;
    if (IsGrid(x, y)) {
        var ij = ScreenToIdx(x, y)
        i0 = ij[0]
        j0 = ij[1]
    }
}

function StopSelect() {
    selectcell = false;
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
}

function StartSelectPos(cr) {
    StopSelect()
    StopSelectPos()

    circrad = cr

    selectpos = true;
    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
}

function ApplySelectPos() {
    // draw shape at (mouse.x mouse.y)
    var xyp = ScreenToPhys(mouse.x, mouse.y)
    var lx = 1, ly = lx*ny/nx
    var dx = lx/nx, dy = ly/ny
    // XXX: implies dx == dy
    var a = circrad * dx;
    var b = circrad * dy;
    var el = {x0: xyp[0] * dx, y0: xyp[1] * dy, a: a, b: b}
    var vof = new Vof(dx, dy, vof_ellipse, el)
    var du = matrix_new(nx, ny)
    vof.grid(nx, ny, /**/ du)
    var i, j
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            u[i][j] = Clip(u[i][j] + du[i][j], 0, 1)
        }
    }
}

function StopSelectPos() {
    selectpos = false;
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
}

function SetSize(nxa, nya) {
    nx = nxa
    ny = nya
    UpdatePar()
    UpdateText(GridToText(u, ny))
    onPaint();
}

function SetButtons() {
    var d = document.getElementById('buttons');
    var bb = d.getElementsByTagName('button');
    var i, bl = bb.length;
    for (i = 0, bl = bb.length; i < bl; ++i) {
        var b = bb[i]
        b.style.width = (0.9 * wx / bl).toString() + "px"
        b.style.height = (0.6 * wx / bl).toString() + "px"
        b.style.fontSize = (wx * 0.05) + "px"
        b.style.padding = "0"
    }
    textarea.style.width = (wx * 0.9) + "px"
    textarea.style.fontSize = (wx * 0.03) + "px"
    ctx.lineWidth = (wx * 0.005)
    wbar = wx* 0.02
}


// Parameters
var nx = 6
var ny = 6
var kbar = 0.2   // bar length relative to screen width
var marx = 0.05  // x-margin relative to screen width
var mary = kbar * 0.5  // y-margin relative to screen width
var base, w, wx, wy, wbar

var canvas = document.getElementById('myCanvas');
var ctx


function UpdatePar() {
    wx = window.innerWidth || document.documentElement.clientWidth || document.body.clientWidth;
    wy = window.innerHeight || document.documentElement.clientHeight || document.body.clientHeight;
    w = (wx * (1.0 - marx * 2)) / nx;   // block size
    if (w * ny > 0.55 * wy) {
        wx = Math.min(wx, wy * 0.6) 
        w = (wx * (1.0 - marx * 2)) / nx;
    }
    base = {x: wx * marx , y: wx * mary}

    u = matrix_zero(nx, ny)

    canvas.width = nx * w + base.x * 2
    canvas.height = ny * w + base.y * 2
    canvas.style.width = canvas.width.toString() + "px"
    canvas.style.height = canvas.height.toString() + "px"

    ctx = canvas.getContext('2d');

    ctx.lineWidth = 5;

    SetButtons();
}

var painting = document.getElementById('paint');
var paint_style = getComputedStyle(painting);
var textarea = document.getElementById("myTextarea");

UpdatePar()


var mouse = {x: 0, y: 0}
var start = {x: 0, y: 0}
var mousepressed = false
var selectcell = false
var selectpos = false

mouse.x = base.x
mouse.y = base.y
start.x = base.x
start.y = base.y



function addListenerMulti(element, eventNames, listener) {
    var events = eventNames.split(' ')
    var i, ie;
    for (i = 0, ie = events.length; i < ie; i++) {
        element.addEventListener(events[i], listener, false);
    }
}

addListenerMulti(canvas, 'mousemove mousestart', function(e) {
    mouse.x = e.pageX - this.offsetLeft;
    mouse.y = e.pageY - this.offsetTop;
});

addListenerMulti(canvas, 'touchmove touchstart', function(e) {
    mouse.x = e.touches[0].pageX - this.offsetLeft;
    mouse.y = e.touches[0].pageY - this.offsetTop;
});

addListenerMulti(canvas, 'mousedown touchstart', function(e) {
    if (selectcell) {
        ApplySelect()
        StopSelect()
    } else if (selectpos) {
        ApplySelectPos()
        StopSelectPos()
    } else {
        mousepressed = true; // start dragging cell
        ustart = matrix_copy(nx, ny, u);
        start.x = mouse.x
        start.y = mouse.y
        canvas.addEventListener('mousemove', onPaint, false);
        canvas.addEventListener('touchmove', onPaint, false);
        onPaint()
    }
});

addListenerMulti(canvas, 'mouseup mouseout touchend touchcancel', function(e) {
    if (mousepressed) {
        canvas.removeEventListener('mousemove', onPaint, false);
        canvas.removeEventListener('touchmove', onPaint, false);
        mousepressed = false;
        onPaint()
    }
});

addListenerMulti(textarea, 'change input keyup onpaste oncut', function(e) {
    UpdateGrid()
    Clear()
    DrawAll()
}, false);

function TwoEllipses() {
    var lx = 1, ly = lx*ny/nx
    var u = matrix_new(nx, ny)
    var dx = lx/nx, dy = ly/ny
    var Param = {}
    var ellipseA = {x0: 0.8, y0: 0.5, a: 0.2, b: 0.4}
    var ellipseB = {x0: 0.2, y0: 0.5, a: 0.2, b: 0.4}
    Param.param = [ellipseA, ellipseB]
    Param.f = [vof_ellipse, vof_ellipse]
    var vof = new Vof(dx, dy, vof_comosite, Param)
    vof.grid(nx, ny, /**/ u)
    return u
}

function OneEllipse() {
    var lx = 1, ly = lx*ny/nx
    var u = matrix_new(nx, ny)
    var dx = lx/nx, dy = ly/ny
    var ellipse = {x0: 0.5, y0: 0.5, a: 0.2, b: 0.4}
    var vof = new Vof(dx, dy, vof_ellipse, ellipse)
    vof.grid(nx, ny, /**/ u)
    return u
}

var i0 = 2, j0 = 3
var circrad = 2
var ends

SetSize(nx, ny)

var u = TwoEllipses()
//var u = OneEllipse()

var ustart = matrix_copy(nx, ny, u)

onPaint()
