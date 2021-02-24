var g_lines;
var g_lines_ptr;
var AddVelocityAngle;
var SetRuntimeConfig;
var GetLines;
var Spawn;
var TogglePause;
var g_tmp_canvas;
var kScale = 1;

function GetExtraConfig() {
  return `
set double rho1 1
set double rho2 0.1

set int nsteps 2

set double cfl 0.5

set double growth_rate 100

set double resist1 1
set double resist2 100

set int visinterp 1

set double cmax 0.6
set double surface_height 0.8

set double visvort 0.01
set double vispot 0
set double vistracer 2
set double visvf 0.7

#set double vistracer 0

set double nucleate_cmax_factor 0.7
set double nucleate_noise 0.001
`
}

function SetExtraConfig(conf) {
  Module.ccall('SetExtraConfig', null, ['string'], [conf]);
}

function SetRuntimeConfig(conf) {
  Module.ccall('SetRuntimeConfig', null, ['string'], [conf]);
}

function SetSigma(sigma) {
  let c = `
set double sigma ${sigma}
`;
  SetRuntimeConfig(c);
  return c;
}

function SetMu(mu) {
  let c = `
set double mu1 ${mu}
set double mu2 ${mu * 0.1}
`;
  SetRuntimeConfig(c);
  return c;
}

function SetRate(r) {
  let c = `
set double reaction_rate_left ${r * 2}
set double reaction_rate_right ${r}
`;
  SetRuntimeConfig(c);
  return c;
}

function SetDiffusion(d) {
  let c = `
set vect tracer_diffusion ${d}
`;
  SetRuntimeConfig(c);
  return c;
}

function SetGravity(flag) {
  let g = flag ? -10 : 0;
  let c = `
set vect gravity 0 ${g}
`;
  SetRuntimeConfig(c);
  return c;
}


function Init(nx) {
  conf = GetExtraConfig()
  let cc = [];
  cc.push(SetSigma(window.range_sigma.value));
  cc.push(SetMu(window.range_mu.value));
  cc.push(SetRate(window.range_rate.value));
  cc.push(SetDiffusion(window.range_diffusion.value));
  cc.push(SetGravity(window.checkbox_gravity.checked));
  conf += cc.join('');
  SetExtraConfig(conf);
  Module.ccall('SetMesh', null, ['number'], [nx])
}


function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  ctx.lineWidth = 3;
  ctx.strokeStyle="#000000";
  ctx.strokeRect(0, 0, canvas.width, canvas.height);

  g_lines = new Uint16Array(Module.HEAPU8.buffer, g_lines_ptr, g_lines_max_size);
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
  g_lines_max_size = 10000;
  g_lines_ptr = Module._malloc(g_lines_max_size * 2);
  TogglePause = Module.cwrap('TogglePause', null, []);
  Spawn = Module.cwrap('Spawn', null, ['number', 'number', 'number']);
  AddVelocityAngle = Module.cwrap('AddVelocityAngle', null, ['number']);
  GetLines = Module.cwrap('GetLines', 'number', ['number', 'number']);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = canvas.width / kScale;
  g_tmp_canvas.height = canvas.height / kScale;

  let keydown = function(ev){
    if (ev.key == ' ') {
      TogglePause();
      ev.preventDefault();
    }
  };
  let keyup = function(ev){
    if (ev.key == ' ') {
      ev.preventDefault();
    }
  };
  let mouseclick = function(ev){
    let x = ev.offsetX / canvas.width;
    let y = 1 - ev.offsetY / canvas.height;
    Spawn(x, y, 0.05);
  };

  window.addEventListener('keydown', keydown, false);
  window.addEventListener('keyup', keyup, false);
  canvas.addEventListener('click', mouseclick, false);
  Init(32);
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
