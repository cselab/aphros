function Clip(a,l,h) {
    if (isFinite(a)) {
        return Math.min(Math.max(a, l), h);
    } 
    return 0.
}

function sqr(a) {
    return a * a;
}

function InitGrid(nx, ny) {
    var u = [];

    for (var j = 0; j < ny; j++) {
        u[j] = [];
        for (var i = 0; i < nx; i++) {
            u[j][i] = 1. / (sqr(i - (nx - 1) * 0.5) * 1.1 +
                sqr(j - (ny - 1) * 0.5) * 1.3 + 0.01);
            if (u[j][i] < 0.3) {
                u[j][i] = 0.;
            }
            u[j][i] = Clip(u[j][i], 0., 1.)
        }
    }

    return u;
}

function CopyGrid(nx, ny, us) {
    var u = [];

    for (var j = 0; j < ny; j++) {
        u[j] = [];
        for (var i = 0; i < nx; i++) {
            u[j][i] = us[j][i];
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
    var k = 2.;
    var wk = w * k;

    var i = Clip(Math.floor((start.x - base.x) / w), 0, nx - 1);
    var j = Clip(Math.floor((start.y - base.y) / w), 0, ny - 1);

    // displacement [px]
    dy = mouse.y - start.y;
    // range of u correction: [-dum,dup]
    dup = 1. - ustart[j][i]
    dum = ustart[j][i]
    // normalize
    dyn = dy / wk

    dyp = Finv(dup) * wk
    dym = Finv(dum) * wk

    dy = Clip(dy, -dyp, dym)

    // correction
    du = -F(Math.abs(dyn)) * Math.sign(dyn);

    u[j][i] = Clip(ustart[j][i] + du, 0., 1.);
    return [dy, dym, dyp]
}

function DrawGrid(u) {
    for (var j = 0 ; j < ny; j++) {
        for (var i = 0 ; i < nx; i++) {
            x = base.x + i * w;
            y = base.y + j * w;
            ctx.strokeRect(x, y, w, w);
            q = u[j][i];
            q = Math.floor((1-q) * 255 + q * 127);
            ctx.fillStyle = RGBToHex(q, q, q);
            ctx.fillRect(x, y, w, w);

            ctx.font = "20px Arial";
            ctx.fillStyle = "#000000";
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            ctx.fillText((u[j][i]).toFixed(2), x + w * 0.5, y + w * 0.5);
        }
    }
}

function DrawLines(ee, u) {
    for (var j = 0 ; j < ny; j++) {
        for (var i = 0 ; i < nx; i++) {
            if (u[i][j] <= 0 || u[i][j] >= 1) {
                continue
            }
            e = ee[i][j]
            xq = nx
            yq = ny
            xa = base.x + (e[1]) * w
            ya = base.y + (e[0]) * w
            xb = base.x + (e[3]) * w
            yb = base.y + (e[2]) * w

            // unit normal
            mx = yb - ya
            my = -(xb - xa)
            mm = Math.sqrt(mx * mx + my * my)
            mx /= mm
            my /= mm
            // shift [px]
            sh = 1
            dx = sh * mx
            dy = sh * my

            ctx.strokeStyle = "#808080"
            ctx.beginPath();
            ctx.moveTo(xa + dx, ya + dy)
            ctx.lineTo(xb + dx, yb + dy)
            ctx.closePath()
            ctx.stroke()

            ctx.strokeStyle = "#000000"
            ctx.beginPath();
            ctx.moveTo(xa - dx, ya - dy)
            ctx.lineTo(xb - dx, yb - dy)
            ctx.closePath()
            ctx.stroke()
        }
    }
}

function DrawString(pp) {
    ctx.beginPath();
    for (var i = 1 ; i < pp.length; ++i) {
        xa = base.x + (pp[i-1][0]) * w
        ya = base.y + (pp[i-1][1]) * w
        xb = base.x + (pp[i][0]) * w
        yb = base.y + (pp[i][1]) * w

        ctx.strokeStyle = "red";
        ctx.beginPath();
        ctx.moveTo(xa, ya)
        ctx.lineTo(xb, yb)
        ctx.closePath();
        ctx.stroke()
    }

    for (var i = 0 ; i < pp.length; ++i) {
        x = base.x + (pp[i][0]) * w
        y = base.y + (pp[i][1]) * w

        ctx.strokeStyle = "black";
        ctx.beginPath();
        ctx.arc(x, y, 3, 0, 2 * Math.PI);
        ctx.closePath();
        ctx.stroke();
    }
    ctx.closePath();
}

function DrawBar(dy, dym, dyp) {
    q = 8
    qh = q / 2
    qq = q / 2
    qqh = qq / 2
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
    for (var j = 0; j < ny; ++j) {
        t += (j == 0 ? "" : "\n")
        for (var i = 0; i < nx; ++i) {
            t += (i == 0 ? "" : " ") + u[j][i]
        }
    }
    return t
}

function TextToGrid(t, nx, ny, /**/ u) {
    t = t.replace(/\s+/g, " ")
    var v = t.split(" ")
    m = 0
    for (var j = 0; j < ny; j++) {
        for (var i = 0; i < nx; i++) {
            u[j][i] = Clip(parseFloat(v[m++]), 0., 1.)
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

function DrawAll() {
    matrix_halo_zero(nx, ny, 2, u)
    ends = matrix_new(nx, ny)
    partstr_vof_ends(nx, ny, u, ends)

    DrawGrid(u)
    DrawLines(ends, u)
    DrawString([[0.3, 0.4], [1, 1], [2, 2], [3, 2.5]])
}

function Clear() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
}

var onPaint = function() {
    Clear()

    if (mousepressed) {
        dd = IncCell()
    }

    DrawAll()

    if (mousepressed) {
        DrawBar(dd[0], dd[1], dd[2])
    }

    UpdateText(GridToText(u, ny))
};


// Parameters
var w = 60;   // block size
var nx = 5
var ny = 5
var base = {x: 50, y: 200}

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

mouse.x = base.x
mouse.y = base.y
start.x = base.x
start.y = base.y

var painting = document.getElementById('paint');
var paint_style = getComputedStyle(painting);

var textarea = document.getElementById("myTextarea");

canvas.addEventListener('mousemove', function(e) {
    mouse.x = e.pageX - this.offsetLeft;
    mouse.y = e.pageY - this.offsetTop;
}, false);

canvas.addEventListener('touchmove', function(e) {
    x = e.touches[0].pageX - this.offsetLeft;
    y = e.touches[0].pageY - this.offsetTop;
    ctx.fillStyle = "green"
    ctx.fillRect(x, y, w, w);
}, false);

canvas.addEventListener('mousedown', function(e) {
    mousepressed = true
}, false);

canvas.addEventListener('mouseup', function(e) {
    mousepressed = false
}, false);

canvas.addEventListener('mouseout', function(e) {
    mousepressed = false
}, false);

canvas.addEventListener('mousedown', function(e) {
    start.x = mouse.x
    start.y = mouse.y

    ustart = CopyGrid(nx, ny, u);

    onPaint()
    canvas.addEventListener('mousemove', onPaint, false);
}, false);

canvas.addEventListener('mouseup', function() {
    onPaint()
    canvas.removeEventListener('mousemove', onPaint, false);
}, false);

canvas.addEventListener('mouseout', function() {
    canvas.removeEventListener('mousemove', onPaint, false);
}, false);

textarea.addEventListener('change', function(e) {
    UpdateGrid()
    Clear()
    DrawAll()
}, false);


u = InitGrid(nx, ny);
ustart = CopyGrid(nx, ny, u);
onPaint();
