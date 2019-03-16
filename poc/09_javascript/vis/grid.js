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
function InitGrid(nx, ny) {
    var u = []
    var i, j

    for (i = 0; i < nx; i++) {
        u[i] = [];
        for (j = 0; j < ny; j++) {
            u[i][j] = 1.0 / (sqr(i - (nx - 1) * 0.5) * 1.1 +
                sqr(j - (ny - 1) * 0.5) * 1.3 + 0.01);
            if (u[i][j] < 0.3) {
                u[i][j] = 0.0;
            }
            u[i][j] = Clip(u[i][j], 0.0, 1.0)
        }
    }

    return u;
}

function CopyGrid(nx, ny, us) {
    var u = [];
    var i, j

    for (i = 0; i < nx; i++) {
        u[i] = [];
        for (j = 0; j < ny; j++) {
            u[i][j] = us[i][j];
        }
    }

    return u;
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
    var k = 2.0;
    var wk = w * k;

    var ij = ScreenToIdx(start.x, start.y)
    var i = ij[0]
    var j = ij[1]

    // displacement [px]
    var dy = mouse.y - start.y;
    // range of u correction: [-dum,dup]
    var dup = 1.0 - ustart[i][j]
    var dum = ustart[i][j]
    // normalize
    var dyn = dy / wk

    var dyp = Finv(dup) * wk
    var dym = Finv(dum) * wk

    dy = Clip(dy, -dyp, dym)

    // correction
    var du = -F(abs(dyn)) * sign(dyn);

    u[i][j] = Clip(ustart[i][j] + du, 0.0, 1.0);
    return [dy, dym, dyp]
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
            var sh = 1
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

function DrawString(pp) {
    ctx.beginPath();
    var i
    for (i = 1 ; i < pp.length; ++i) {
        var xya = PhysToScreen(pp[i - 1][0], pp[i - 1][1])
        var xyb = PhysToScreen(pp[i][0], pp[i][1])
        var xa = xya[0], ya = xya[1]
        var xb = xyb[0], yb = xyb[1]

        ctx.strokeStyle = "red";
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

function DrawBar(dy, dym, dyp) {
    var q = 8
    var qh = q / 2
    var qq = q / 2
    var qqh = qq / 2
    ctx.fillStyle = '#a0a0a0';
    ctx.fillRect(start.x-qh, start.y +dym, q, -dyp-dym-qq);
    ctx.fillStyle = 'blue';
    ctx.fillRect(start.x-qqh, start.y-qqh, qq, dy);
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
            t += (i == 0 ? "" : " ") + u[i][ny - j - 1].toFixed(4)
        }
    }
    return t
}

function TextToGrid(t, nx, ny, /**/ u) {
    t = t.replace(/\s+/g, " ")
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
    t = textarea.value
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

function DrawAll() {
    var AX = 0, AY = 1, BX = 2, BY = 3

    matrix_halo_zero(nx, ny, 2, u)
    var ends = matrix_new(nx, ny)
    partstr_vof_ends(nx, ny, u, ends)
    DrawGrid(u)
    DrawSelectOld(i0, j0)
    DrawLines(ends, u)
    if (ends[i0][j0] !== undefined) {
        var end = partstr_cell_ends(nx, ny, i0, j0, ends)
        var ne = end.length
        var nh = 4
        var hp = 4.0 / (2.0*nh)
        var eta = 0.5
        var t = 0.0
        var n = 2*nh + 1
        var partstr = new Partstr(nh, hp, eta)
        var e = end[0]
        var p = [(e[AX] + e[BX]) * 0.5, (e[AY] + e[BY]) * 0.5]
        var dx = e[BX] - e[AX]
        var dy = e[BY] - e[AY]
        var a = Math.atan2(dy, dx)
        partstr.start(ne, end, a, t, p)
        for (var k = 0; k < 100; k++)
            partstr.step()
        DrawString(partstr.xx)
    }
}

function Clear() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
}

var onPaint = function() {
    Clear()


    if (mousepressed && selectcell) {
        var ij = ScreenToIdx(start.x, start.y)
        i0 = ij[0]
        j0 = ij[1]
        StopSelect()
        DrawAll()
        return
    }

    if (mousepressed) {
        var dd = IncCell()
    }

    DrawAll()

    if (mousepressed) {
        DrawBar(dd[0], dd[1], dd[2])
    }

    if (selectcell) {
        var ij = ScreenToIdx(mouse.x, mouse.y)
        DrawSelectCurrent(ij[0], ij[1], 1)
    }

    UpdateText(GridToText(u, ny))
};

function StartSelect() {
    selectcell = true;
    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
}

function StopSelect() {
    selectcell = false;
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
}


// Parameters
var nx = 5
var ny = 5
var marx = 0.1  // x-margin relative to screen width
var maryb = 1  // y-margin relative to block size
var wx = window.innerWidth || document.documentElement.clientWidth || document.body.clientWidth;
wx = Math.min(wx, 1000)
var w = (wx * (1.0 - marx * 2)) / nx;   // block size
var base = {x: wx * marx , y: w * maryb}

var canvas = document.getElementById('myCanvas');

canvas.width = nx * w + base.x * 2
canvas.height = ny * w + base.y * 2
canvas.style.width = canvas.width.toString() + "px"
canvas.style.height = canvas.height.toString() + "px"


var ctx = canvas.getContext('2d');

ctx.lineWidth = 3;
ctx.lineJoin = 'round'
ctx.lineCap = 'round'
ctx.strokeStyle = '#505050'

var mouse = {x: 0, y: 0}
var start = {x: 0, y: 0}
var mousepressed = false
var selectcell = false

mouse.x = base.x
mouse.y = base.y
start.x = base.x
start.y = base.y

var painting = document.getElementById('paint');
var paint_style = getComputedStyle(painting);

var textarea = document.getElementById("myTextarea");


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
    mousepressed = true
});

addListenerMulti(canvas, 'mouseup mouseout touchend touchcancel', function(e) {
    mousepressed = false
});

addListenerMulti(canvas, 'mousedown touchstart', function(e) {
    start.x = mouse.x
    start.y = mouse.y

    ustart = CopyGrid(nx, ny, u);

    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
});

addListenerMulti(canvas, 'mouseup mouseout touchend touchcancel', function(e) {
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
});

textarea.addEventListener('change', function(e) {
    UpdateGrid()
    Clear()
    DrawAll()
}, false);


var u = InitGrid(nx, ny);
var ustart = CopyGrid(nx, ny, u);
var i0 = 1, j0 = 3;
onPaint();
