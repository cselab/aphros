var g_lines;
var g_lines_ptr;
var SetRuntimeConfig;
var Spawn;
var g_tmp_canvas;
var kScale = 1;

function GetExtraConfig() {
  return `
set double cfl 0.9
set double cflvis 0.125
set double cflsurf 2

set double tmax 100
set double vispressure 1
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
set double mu2 ${mu * 0.01}
`;
  SetRuntimeConfig(c);
  return c;
}

function SetGravity(g) {
  let c = `
set vect gravity 0 ${-g}
`;
  SetRuntimeConfig(c);
  return c;
}

function TogglePause() {
  let s = Module.ccall('TogglePause', 'number', []);
  let button = document.getElementById('button_pause');
  if (button) {
    if (s) {
      button.className = "button pressed";
    } else {
      button.className = "button";
    }
  }
}

function ResetButtons() {
  [16, 32, 64, 128].forEach(nx => {
    let button = document.getElementById('button_' + nx);
    if (button) {
      button.className = "button";
    }
  });
}

function Init(nx) {
  conf = GetExtraConfig()
  let cc = [];
  if (window.checkbox_coal.checked) {
    cc.push("set int coal 1\n");
  } else {
    cc.push("set int coal 0\n");
  }
  cc.push(SetSigma(window.range_sigma.value));
  cc.push(SetMu(window.range_mu.value));
  cc.push(SetGravity(window.range_gravity.value));
  ResetButtons();
  let button = document.getElementById('button_' + nx);
  if (button) {
    button.className = "button pressed";
  }
  conf += cc.join('');
  SetExtraConfig(conf);
  Module.ccall('SetMesh', null, ['number'], [nx])
  //Spawn(0.5, 0.5, 0.2);
}


function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  ctx.lineWidth = 3;
  ctx.strokeStyle="#000000";
  ctx.strokeRect(0, 0, canvas.width, canvas.height);
}

function PostRun() {
  window.checkbox_coal.checked = false;
  g_lines_max_size = 10000;
  g_lines_ptr = Module._malloc(g_lines_max_size * 2);
  Spawn = Module.cwrap('Spawn', null, ['number', 'number', 'number']);

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
    Spawn(x, y, 0.1);
  };

  window.addEventListener('keydown', keydown, false);
  window.addEventListener('keyup', keyup, false);
  canvas.addEventListener('click', mouseclick, false);
  Init(16);
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
