var g_lines;
var g_lines_ptr;
var SetRuntimeConfig;
var GetLines;
var Spawn;
var g_tmp_canvas;
var kScale = 1;

function GetExtraConfig() {
  let res = `
set string Cred 1 0.12 0.35
set string Cgreen 0 0.8 0.42
set string Cblue 0 0.6 0.87
set string Cwhite 1 1 1
set string Cblack 0 0 0

set double cflvis 1000
set double cflst 1

set string visual "
vf {
set vect values 0 1
set vect colors $Cwhite $Cgreen
set vect opacities 1 1
}
"
`
  return res;
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
  let input_conf = document.getElementById('input_conf');
  if (input_conf) {
    cc.push(input_conf.value);
  }
  ResetButtons();
  let button = document.getElementById('button_' + nx);
  if (button) {
    button.className = "button pressed";
  }
  conf += cc.join('');
  SetExtraConfig(conf);
  Module.ccall('SetMesh', null, ['number'], [nx])
}


function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, 0, 0, canvas.width, canvas.height);

  // Draw bounding box
  ctx.lineWidth = 3;
  ctx.strokeStyle="#000000";
  ctx.strokeRect(0, 0, canvas.width, canvas.height);

  g_lines = new Uint16Array(Module.HEAPU8.buffer, g_lines_ptr, g_lines_max_size);

  // Draw interface lines
  {
    let size = GetLines(0, g_lines.byteOffset, g_lines.length);
    ctx.lineWidth = 3;
    ctx.strokeStyle = "black";
    for (let i = 0; i + 3 < size; i += 4) {
      ctx.beginPath();
      ctx.moveTo(g_lines[i], g_lines[i + 1]);
      ctx.lineTo(g_lines[i + 2], g_lines[i + 3]);
      ctx.stroke();
    }
  }

  // Draw embed lines
  {
    let size = GetLines(1, g_lines.byteOffset, g_lines.length);
    ctx.lineWidth = 3;
    ctx.strokeStyle = "black";
    for (let i = 0; i + 3 < size; i += 4) {
      ctx.beginPath();
      ctx.moveTo(g_lines[i], g_lines[i + 1]);
      ctx.lineTo(g_lines[i + 2], g_lines[i + 3]);
      ctx.stroke();
    }
  }
}

function PostRun() {
  g_lines_max_size = 10000;
  g_lines_ptr = Module._malloc(g_lines_max_size * 2);
  Spawn = Module.cwrap('Spawn', null, ['number', 'number', 'number']);
  GetLines = Module.cwrap('GetLines', 'number', ['number', 'number', 'number']);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = canvas.width / kScale;
  g_tmp_canvas.height = canvas.height / kScale;

  let keydown = function(ev){
    if (ev.key == ' ') {
      TogglePause();
      //ev.preventDefault();
    }
  };
  let keyup = function(ev){
    if (ev.key == ' ') {
      //ev.preventDefault();
    }
  };
  let mouseclick = function(ev){
    let x = ev.offsetX / canvas.width;
    let y = 1 - ev.offsetY / canvas.height;
    Spawn(x, y, 0.1);
  };

  canvas.addEventListener('click', mouseclick, false);
  window.addEventListener('keydown', keydown, false);
  window.addEventListener('keyup', keyup, false);
  [
    window.input_conf,
  ].forEach(b => {
    b.addEventListener('keydown', function(ev){ev.stopPropagation();}, false);
    b.addEventListener('keyup', function(ev){ev.stopPropagation();}, false);
  });
  [
    window.button_pause, window.button_16, window.button_32,
    window.button_64, window.button_128, 
  ].forEach(b => {
    b.addEventListener('keydown', function(ev){
      if (ev.key == ' ') {
        ev.preventDefault();
      }
    }, false);
    b.addEventListener('keyup', function(ev){
      if (ev.key == ' ') {
        ev.preventDefault();
      }
    }, false);
  });
  Init(64);
}

var Module = {
  preRun: [],
  postRun: [PostRun],
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
  printErr: (function(text) {
    var element = document.getElementById('outputerr');
    if (element) element.value = '';
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      console.error(text);
      if (element) {
        element.value += text + "\n";
        element.scrollTop = element.scrollHeight;
      }
    };
  })(),
  canvas: (function() { return document.getElementById('canvas'); })(),
  setStatus: function(text) {},
};
