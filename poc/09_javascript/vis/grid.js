var canvas = document.getElementById('myCanvas');
var ctx = canvas.getContext('2d');

var painting = document.getElementById('paint');
var paint_style = getComputedStyle(painting);
canvas.width = parseInt(paint_style.getPropertyValue('width'));
canvas.height = parseInt(paint_style.getPropertyValue('height'));

var mouse = {x: 0, y: 0};
var base = {x: 100, y: 100}

mouse.x = base.x
mouse.y = base.y
canvas.addEventListener('mousemove', function(e) {
  mouse.x = e.pageX - this.offsetLeft;
  mouse.y = e.pageY - this.offsetTop;
}, false);

function Clip(a,l,h) {
    return Math.min(Math.max(a, l), h);
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

    canvas.addEventListener('mousemove', onPaint, false);
}, false);

canvas.addEventListener('mouseup', function() {
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
    return Math.pow(a, 2);
}


var onPaint = function() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    var w = 60;

    var i = Clip(Math.floor((start.x - base.x) / w), 0, nx - 1);
    var j = Clip(Math.floor((start.y - base.y) / w), 0, ny - 1);

    dy = start.y - mouse.y;
    dy = dy / w;
    dy = F(Math.abs(dy)) * Math.sign(dy) / F(2.);
    u[j][i] = Clip(ustart[j][i] + dy, 0., 1.);

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
    ctx.fillStyle = 'red';
    ctx.fillRect(mouse.x, mouse.y, 10, 10);
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
