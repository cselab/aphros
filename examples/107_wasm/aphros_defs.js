var g_lines;
var g_lines_ptr;
var SetRuntimeConfig;
var GetLines;
var Spawn;
var g_tmp_canvas;
var kScale = 1;
var input_conf = document.getElementById('input_conf');
var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');

function EncodeSafe(base) {
  return base
          .replace(/\+/g, '-')
          .replace(/\//g, '_');
}

function DecodeSafe(base) {
  return base
          .replace(/-/g, '+')
          .replace(/_/g, '/');
}

function GetExtraConfig() {
  let res = `
set string Cred 1 0.12 0.35
set string Cgreen 0 0.8 0.42
set string Cblue 0 0.6 0.87
set string Cwhite 1 1 1
set string Cblack 0 0 0
set string Cgray 0.5 0.5 0.5

set string V_volume_fraction_green "
volume fraction {
set vect values 0 1
set vect colors $Cgreen $Cgreen
set vect opacities 0 1
}
"

set string V_embed_gray "
volume fraction {
set vect values 0 1
set vect colors $Cgray $Cgray
set vect opacities 1 0
}
"

`
  return res;
}

function Compress(text) {
  return EncodeSafe(LZString.compressToBase64(text));
}

function Decompress(compressed) {
  return LZString.decompressFromBase64(DecodeSafe(compressed));
}

function SetExtraConfig(config) {
  Module.ccall('SetExtraConfig', null, ['string'], [config]);
}

function ApplyConfig() {
  UpdateFullUrl();
  SetRuntimeConfig(input_conf.value);
}

function SetRuntimeConfig(config) {
  Module.ccall('SetRuntimeConfig', null, ['string'], [config]);
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
  ClearOutput();
  config = GetExtraConfig()
  let cc = [];
  UpdateFullUrl();
  cc.push(input_conf.value);
  ResetButtons();
  let button = document.getElementById('button_' + nx);
  if (button) {
    button.className = "button pressed";
  }
  config += cc.join('');
  SetExtraConfig(config);
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

function ClearOutput() {
  output.value = '';
  outputerr.value = '';
}
function Print(text) {
  //console.log(text);
  if (output) {
    output.value += text + "\n";
    output.scrollTop = output.scrollHeight;
  }
}
function PrintError(text) {
  console.error(text);
  if (outputerr) {
    outputerr.value += text + "\n";
    outputerr.scrollTop = outputerr.scrollHeight;
  }
}
function GetFullUrl(config) {
  let url = new URL(location.href);
  url.search = "?config=" + Compress(config);
  return url.toString();
}
function ClearUrl() {
  let url = new URL(location.href);
  url.search = "";
  history.pushState(null, null, url.toString());
}
function UpdateFullUrl() {
  let fullurl = GetFullUrl(input_conf.value);
  history.pushState(null, null, fullurl);
  let element = document.getElementById('a_fullurl');
  if (element) {
    element.href = fullurl;
  }
}
function Request(url, action) {
  let req = new XMLHttpRequest();
  req.onload = function() {
    if (req.status == 200) {
      action(req.response)
    } else {
      PrintError(`Error ${req.status}: ${req.statusText}`);
    }
  };
  req.open('GET', url, true);
  req.send();
}
function GetShortUrl(fullurl, action) {
  return Request("https://tinyurl.com/api-create.php?url=" + fullurl, action);
}
function UpdateShortUrl() {
  let fullurl = GetFullUrl(input_conf.value);
  GetShortUrl(fullurl, function(response) {
    document.getElementById('text_shorturl').value = response;
  });
}

function PostRun() {
  let url = new URL(location.href);
  let compressed = url.searchParams.get('config');
  if (compressed) {
    input_conf.value = Decompress(compressed);
  }
  UpdateFullUrl();
  ClearUrl();

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
    input_conf,
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
    ClearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      Print(text);
    };
  })(),
  printErr: (function(text) {
    ClearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      PrintError(text);
    };
  })(),
  canvas: (function() { return document.getElementById('canvas'); })(),
  setStatus: function(text) {},
};
