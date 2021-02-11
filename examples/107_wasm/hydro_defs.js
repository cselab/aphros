var g_lines;
var TogglePause;
var Spawn;
var AddVelocityAngle;

function SetExtraConfig() {
  Module.ccall('SetExtraConfig', '', ['string'], [
          `
set double rho1 1
set double mu1 0.0001
set double rho2 10
set double mu2 0.001
set vect gravity 0 -5
set double hypre_symm_tol 1e-1
set int hypre_symm_maxiter 20
set int hypre_symm_miniter 1
set double visvel 0
set double visvf 0.5
set int sharpen 1
set double dtmax 0.005
set double sigma 4
set int nsteps 3

set double cflsurf 2
`
        ]);
}

function SetMesh(nx) {
  SetExtraConfig();
  Module.ccall('SetMesh', null, ['number'], [nx])
  Spawn(0.5, 0.5, 0.2);
}


function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.lineWidth = 3;
  ctx.strokeStyle="#000000";
  ctx.strokeRect(0, 0, canvas.width, canvas.height);

  var GetLines = Module.cwrap('GetLines', 'number', ['number', 'number']);
  let size = GetLines(g_lines.byteOffset, g_lines.length);
  ctx.lineWidth = 3;
  ctx.strokeStyle = "black";
  for (let i = 0; i + 3 < size; i += 4) {
    ctx.beginPath();
    ctx.moveTo(g_lines[i], g_lines[i + 1]);
    ctx.lineTo(g_lines[i + 2], g_lines[i + 3]);
    ctx.stroke();
  }
}

function PostRun() {
  let max_size = 10000;
  g_lines = new Uint16Array(
    Module.HEAPU8.buffer, Module._malloc(max_size * 2), max_size);
  TogglePause = Module.cwrap('TogglePause', null, []);
  Spawn = Module.cwrap('Spawn', null, ['number', 'number', 'number']);
  AddVelocityAngle = Module.cwrap('AddVelocityAngle', null, ['number']);
  let keypress = function(ev){
    if (ev.key == ' ') {
      TogglePause();
    }
    return false;
  };
  let mouseclick = function(ev){
    let x = ev.offsetX / 500;
    let y = 1 - ev.offsetY / 500;
    Spawn(x, y, 0.1);
    return false;
  };
  window.addEventListener('keypress', keypress, false);
  let canvas = Module['canvas'];
  canvas.addEventListener('click', mouseclick, false);
  SetMesh(32);
}

var Module = {
  preRun: [],
  postRun: [function (){ PostRun();}],
  print: (function() {
    var element = document.getElementById('output');
    if (element) element.value = '';
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      console.log(text);
      if (element) {
        element.value += text + "\n";
        element.scrollTop = element.scrollHeight;
      }
    };
  })(),
  printErr: function(text) {
    if (arguments.length > 1) {
      text = Array.prototype.slice.call(arguments).join(' ');
    }
    console.error(text);
  },
  canvas: (function() { return document.getElementById('canvas'); })(),
  setStatus: function(text) {},
};
