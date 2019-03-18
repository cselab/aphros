// naming conventions
// x,y: screen coordinates
// xp,yp: physical space
// i,j: indices on the grid


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

function ValidScreen(x, y) {
    var i, j;
    i = Math.floor((x - base.x) / w);
    j = ny - Math.floor((y - base.y) / w) - 1;
    return i >= 0 && i < nx && j >= 0 && j < ny;
}

function ValidIdx(i, j) {
    return i >= 0 && i < nx && j >= 0 && j < ny
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

function IncCell() {
    var wk = wx * barl_

    var ij = ScreenToIdx(start.x, start.y)
    var i = ij[0]
    var j = ij[1]

    var dirx = (abs(mouse.x - start.x) >= abs(mouse.y - start.y));
    // range of u correction: [dum,dup]
    var dup = 1.0 - ustart[i][j]
    var dum = -ustart[i][j]
    // displacement [px]
    var dy = (dirx ? mouse.x - start.x : -(mouse.y - start.y))
    // normalize
    var dyn = dy / wk

    var dyp = dup * wk
    var dym = dum * wk

    dy = Clip(dy, dym, dyp)

    // correction
    var du = dy / wk

    u[i][j] = Clip(ustart[i][j] + du, 0.0, 1.0);
    return [dy, dym, dyp, dirx]
}


function DrawCellFill(q, x, y) {
    ctx.strokeRect(x, y, w, w);
    q = Math.floor((1-q) * 255 + q * 127);
    ctx.fillStyle = RGBToHex(q, q, q);
    ctx.fillRect(x, y, w, w);
}

function DrawCellText(q, x, y) {
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
            DrawCellFill(q, xy[0], xy[1])
            /* DrawCellText(q, xy[0], xy[1]) */
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
}

function DrawBar(dy, dym, dyp, dirx) {
    var x = start.x, y = start.y;
    var q = barw_ * wx
    var qh = q / 2
    var qq = q / 2
    var qqh = qq / 2
    ctx.fillStyle = '#a0a0a0';
    if (dirx) {
        ctx.fillRect(x+dym-qqh, y-qh, dyp-dym+qh, q);
    } else {
        ctx.fillRect(x-qh, y-dyp-qqh, q, dyp-dym+qh);
    }
    ctx.fillStyle = 'blue';
    if (dirx) {
        ctx.fillRect(x, y-qqh, dy, qq);
    } else {
        ctx.fillRect(x-qqh, y, qq, -dy);
    }
    ctx.fillStyle = 'green';
    ctx.fillRect(x-qh, y-qh, q, q);
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

function SplitText(t) {
    var ll = []
    var ss = t.split(/[\r\n]+/g);
    var ns = ss.length;
    var j;
    for (j = 0; j < ns; ++j) {
        var s = ss[j].replace(/\s+/g, ' ')
        s = s.replace(/^\s+|\s+$/g, '')
        ll[j] = s.split(' ')
    }
    return ll
}

function TextToSize(t) {
    var ll = SplitText(t)
    var ny = ll.length;
    var nx = 0;
    if (ny > 0) {
        var j;
        for (j = 0; j < ny; ++j) {
            nx = Math.max(nx, ll[j].length)
        }
    }
    return [nx, ny]
}

function TextToGrid(t, nx, ny, /**/ u) {
    var ll = SplitText(t)
    var i, j
    for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
            u[i][ny - j - 1] = Clip(parseFloat(ll[j][i]), 0.0, 1.0)
        }
    }
    return u
}

function UpdateText0(t) {
    textarea.value = t
}

function UpdateText() {
    UpdateText0(GridToText(u, ny))
}

function UpdateGrid() {
    var t = textarea.value
    var nxy = TextToSize(t)
    if (nxy[0] != nx || nxy[1] != ny) {
        SetSize0(nxy[0], nxy[1])
    }
    TextToGrid(t, nx, ny, u)
}

function DrawSelect(i, j, next) {
    var xy = IdxToScreen(i, j)
    var x = xy[0] + w/2, y = xy[1] + w/2;
    var q = w / 10
    ctx.fillStyle = (next ? cl_next : cl_curr)
    ctx.fillRect(x - q, y - q, 2 * q, 2 * q);
}

function DrawSelectCurr(i, j) {
    DrawSelect(i, j, false)
}

function DrawSelectNext(i, j) {
    DrawSelect(i, j, true)
}

function DrawSelectShape(xp, yp) {
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
    if (ValidIdx(i0, j0) && ends[i0][j0] !== undefined) {
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
    DrawSelectCurr(i0, j0)
    DrawLines(ends, u)
    DrawStringRun(i0, j0, cl_curr)
}

function Clear() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
}

var onPaint = function() {
    Clear()

    DrawAll()

    if (selectcell) {
        if (ValidScreen(mouse.x, mouse.y)) {
            var ij = ScreenToIdx(mouse.x, mouse.y)
            DrawSelectNext(ij[0], ij[1])
            DrawStringRun(ij[0], ij[1], cl_next)
        }
    } else if (selectshape) {
        var xyp = ScreenToPhys(mouse.x, mouse.y)
        DrawSelectShape(xyp[0], xyp[1])
    } else if (mousepressed) { // dragging cell
        var dd = IncCell()
        DrawBar(dd[0], dd[1], dd[2], dd[3])
        UpdateText()
    }
}

function StartSelect() {
    StopSelect()
    StopSelectShape()

    selectcell = true;
    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
}

function ApplySelect() {
    var x = mouse.x, y = mouse.y;
    if (ValidScreen(x, y)) {
        var ij = ScreenToIdx(x, y)
        SetCell(ij[0], ij[1])
    }
}

function SetCell(i0a, j0a) {
    i0 = i0a
    j0 = j0a
}

function StopSelect() {
    selectcell = false;
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
}

// cr: cicle radius in cells
function StartSelectShape(cr) {
    StopSelect()
    StopSelectShape()

    circrad = cr

    selectshape = true;
    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
    canvas.addEventListener('touchmove', onPaint, false);
}

// xp,yp: center 
// rp: radius
// physical coordinates
function PutCircle(xp, yp, rp) {
    var lx = 1
    var dx = lx / nx
    var dy = dx
    var el = {x0: xp * dx, y0: yp * dy, a: rp * dx, b: rp * dy}
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

function ApplySelectShape() {
    // draw shape at (mouse.x mouse.y)
    var xyp = ScreenToPhys(mouse.x, mouse.y)
    PutCircle(xyp[0], xyp[1], circrad)
    UpdateText()
}

function StopSelectShape() {
    selectshape = false;
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
    canvas.removeEventListener('touchmove', onPaint, false);
}

function SetSize0(nxa, nya) {
    nx = nxa
    ny = nya
    InitFields()
    UpdateViewport()
}

function SetSize(nxa, nya) {
    SetSize0(nxa, nya)
    UpdateText()
    onPaint()
}


// initial
var nx = 7, ny = 7      // grid size
var i0 = 4, j0 = 4      // select cell

// sizes relative to wx
var barl_ = 0.2         // bar length
var barw_ = 0.02        // bar length
var marx_ = barl_ * 0.5 // x-margin
var mary_ = barl_ * 0.5 // y-margin
var hb_ = 0.1           // button height
var wab_ = 0.9          // button row width
var ht_ = 0.3           // text area height
var wt_ = 0.9           // texte area width
var hfb_ = 0.05          // button font height
var hf_ = 0.025          // font height
var wl_ = 0.005          // line width
var cl_curr = "red"
var cl_next = "green"

// global
var u, ustart, ends
var wx                  // screen width
var w                   // cell width
var base                // left top corner of grid
var circrad             // circle radius
var ctx                 // canvas context


function GetWindowWidth(nx, ny) {
    var wx, wy, w
    wx = window.innerWidth || document.documentElement.clientWidth 
        || document.body.clientWidth;
    wy = window.innerHeight || document.documentElement.clientHeight 
        || document.body.clientHeight;
    w = wx * (1 - marx_ * 2) / nx
    w = Math.min(w, 
        wy / (ny + nx * (mary_ * 2 + hb_ + ht_) / (1 - marx_ * 2)))

    wx = w * nx / (1 - marx_ * 2)
    wy = w * ny + wx * (mary_ * 2 + hb_ + ht_)

    return [wx, wy, w]
}

function InitFields() {
    u = matrix_zero(nx, ny)
    ustart = matrix_copy(nx, ny, u)
    ends = matrix_new(nx, ny)
}

function UpdateViewport() {
    var nm = Math.max(nx, ny)
    var ww = GetWindowWidth(nm, nm)

    wx = ww[0]
    w = ww[2]

    base = {x: wx * marx_ , y: wx * mary_}

    canvas.width = w * nx + base.x * 2
    canvas.height = w * ny + base.y * 2
    canvas.style.width = canvas.width.toString() + "px"
    canvas.style.height = canvas.height.toString() + "px"

    ctx = canvas.getContext('2d');

    var d = document.getElementById('buttons');
    d.style.fontSize = (hfb_ * wx) + "px"
    var bb = d.getElementsByTagName('button');
    var i, bl = bb.length;
    for (i = 0, bl = bb.length; i < bl; ++i) {
        var b = bb[i]
        b.style.width = (wab_ * wx / bl).toString() + "px"
        b.style.height = (hb_ * wx).toString() + "px"
    }
    document.getElementById('selectcell').style.color = cl_curr;

    textarea.style.width = (wt_ * wx) + "px"
    var k_ = 0.9 // XXX: factor to avoid vertical scrolling
    textarea.style.height = (ht_ * wx * k_) + "px"
    textarea.style.fontSize = (hf_ * wx) + "px"
    ctx.lineWidth = (wl_ * wx)

    mouse.x = base.x
    mouse.y = base.y
    start.x = base.x
    start.y = base.y
}

function addListenerMulti(element, eventNames, listener) {
    var events = eventNames.split(' ')
    var i, ie;
    for (i = 0, ie = events.length; i < ie; i++) {
        element.addEventListener(events[i], listener, false);
    }
}

function SetEvents() {
    addListenerMulti(canvas, 'mousemove mousedown', function(e) {
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
        } else if (selectshape) {
            ApplySelectShape()
            StopSelectShape()
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

    addListenerMulti(canvas, 'mouseup mouseout touchend touchcancel',
            function(e) {
        if (mousepressed) {
            canvas.removeEventListener('mousemove', onPaint, false);
            canvas.removeEventListener('touchmove', onPaint, false);
            mousepressed = false;
            onPaint()
        }
    });

    addListenerMulti(textarea, 'change input keyup onpaste oncut',
            function(e) {
        UpdateGrid()
        Clear()
        DrawAll()
    }, false);

    addListenerMulti(window, 'resize',
            function(e) {
        UpdateViewport()
        Clear()
        DrawAll()
    }, false);
}

var canvas = document.getElementById('myCanvas');
var textarea = document.getElementById("myTextarea");
var mouse = {x: 0, y: 0}
var start = {x: 0, y: 0}
var mousepressed = false
var selectcell = false
var selectshape = false

SetSize(nx, ny)
PutCircle(nx / 2 - 1.2, ny / 2 + 1.3, 2)
PutCircle(nx / 2 + 2.1, ny / 2 - 1.8, 1)
SetEvents()
onPaint()
