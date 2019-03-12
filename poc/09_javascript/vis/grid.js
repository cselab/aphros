
var canvas = document.getElementById('myCanvas');
var ctx = canvas.getContext('2d');

var painting = document.getElementById('paint');
var paint_style = getComputedStyle(painting);
canvas.width = parseInt(paint_style.getPropertyValue('width'));
canvas.height = parseInt(paint_style.getPropertyValue('height'));

var mouse = {x: 0, y: 0};
var base = {x: 100, y: 100}
var mousepressed = false
var w = 60;   // block size

var textarea = document.getElementById("myTextarea");

mouse.x = base.x
mouse.y = base.y

canvas.addEventListener('mousemove', function(e) {
    mouse.x = e.pageX - this.offsetLeft;
    mouse.y = e.pageY - this.offsetTop;
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

textarea.addEventListener('change', function(e) {
    UpdateGrid()
    DrawGrid(u)
}, false);


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
            u[j][i] = 1. / (sqr(i - (nx - 1) * 0.5) +
                sqr(j - (ny - 1) * 0.5));
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

function DrawLines(ee) {
    for (var j = 0 ; j < ny; j++) {
        for (var i = 0 ; i < nx; i++) {
            e = ee[i][j]
            xq = nx
            yq = ny
            xa = base.x + (e[1]) * w
            ya = base.y + (e[0]) * w
            xb = base.x + (e[3]) * w
            yb = base.y + (e[2]) * w

            ctx.strokeStyle = "#000000";
            ctx.beginPath();
            ctx.moveTo(xa, ya)
            ctx.lineTo(xb, yb)
            ctx.closePath();
            ctx.stroke()
        }
    }
}

function DrawBar() {
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
        t += (j == 0 ? "" : "\n") + u[j].join(" ")
    }
    return t
}

function TextToGrid(t, nx, ny, /**/ u) {
    var v = t.split("\n").join(" ").split(" ")
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

var onPaint = function() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    if (mousepressed) {
        IncCell()
    }

    matrix_halo_zero(nx, ny, 2, u)
    ends = matrix_new(nx, ny)
    partstr_vof_ends(nx, ny, u, ends)

    DrawGrid(u)
    DrawLines(ends)

    if (mousepressed) {
        DrawBar()
    }

    UpdateText(GridToText(u, ny))
};


ctx.lineWidth = 3;
ctx.lineJoin = 'round';
ctx.lineCap = 'round';
ctx.strokeStyle = '#505050';
var nx = 5;
var ny = 5;
u = InitGrid(nx, ny);
ustart = CopyGrid(nx, ny, u);
var start = {x: 0, y: 0};
onPaint();
