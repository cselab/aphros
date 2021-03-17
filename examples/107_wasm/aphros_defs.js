var g_lines;
var g_lines_ptr;
var SetRuntimeConfig;
var GetConfigString;
var GetConfigDouble;
var GetLines;
var Spawn;
var g_tmp_canvas;
var kScale = 1;
var input_conf = document.getElementById('input_conf');
var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');
var g_sliders = {};
var g_sliders_string = "";

function GetInputConfig(fromslider=false) {
  let config = input_conf.value;
  if (fromslider) {
    let m = config.match(/^FROMSLIDER.*/gm);
    if (m) {
      config = m.join('\n');
    } else {
      config = "";
    }
  }
  config = config.replace(/^FROMSLIDER/gm, '');
  return config;
}

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

// colors from
// https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/tree/master/Colours
function GetExtraConfig() {
  let res = `
set string Cred 1 0.12 0.35
set string Cgreen 0 0.8 0.42
set string Cblue 0 0.6 0.87
set string Cpurple 0.686 0.345 0.729
set string Cyellow 1 0.776 0.118
set string Corange 0.949 0.522 0.133
set string Cwhite 1 1 1
set string Cblack 0 0 0
set string Cgray 0.627 0.694 0.729

set string sliders
set double spawn_r 0.05

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

function ApplyConfig() {
  UpdateFullUrl();
  let config = GetInputConfig(false);
  SetRuntimeConfig(config);
  Print(`applied config of ${config.length} characters`);
  UpdateSliders();
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
  cc.push(GetInputConfig(false));
  ResetButtons();
  let button = document.getElementById('button_' + nx);
  if (button) {
    button.className = "button pressed";
  }
  config += cc.join('');
  SetExtraConfig(config);
  Module.ccall('SetMesh', null, ['number'], [nx])
  ApplyConfig();
}

function UpdateSlider(elem) {
  let slider = g_sliders[elem.id];
  let current = document.getElementById(slider.id_current);
  slider.value = elem.value
  current.innerHTML = `${slider.name} (${slider.variable}=${slider.value})`;
  SetRuntimeConfig(`set double ${slider.variable} ${slider.value}`);
  SetRuntimeConfig(GetInputConfig(true));
}

function UpdateSliders() {
  let div_sliders = document.getElementById('div_sliders');
  if (!div_sliders) {
    return;
  }
  let sliders_string = GetConfigString("sliders");
  if (g_sliders_string != sliders_string) {
    g_sliders_string = sliders_string;
    g_sliders = {};
    div_sliders.innerHTML = "";
    sliders_string.split('\n').forEach(str => {
      str = str.trim();
      if (str.length == 0) {
        return;
      }
      let desc = str.split(/ +/);
      let i = 0;
      let variable = desc.length > i ? desc[i++] : desc;
      let min = desc.length > i ? parseFloat(desc[i++]) : 0;
      let max = desc.length > i ? parseFloat(desc[i++]) : 1;
      let step = (max - min) * 0.01;
      let name = desc.length > i ? desc.slice(i++).join(' ') : variable;
      let value = GetConfigDouble(variable);
      if (isNaN(value)) {
        value = (min + max) / 2;
      }
      let id = "slider_" + variable;
      let id_current = id + "_current";
      g_sliders[id] = {
          min:min, max:max, step:step, variable:variable,
          name:name, value:value, id:id, id_current:id_current,
      };
      div_sliders.innerHTML += `
      <div class="emscripten" style="padding: 0px;">
        <label class="label"><input type="range" class="slider" min="${min}" max="${max}" value="${value}" step="${step}" id="${id}"
          onchange="UpdateSlider(this);">
          <span id="${id_current}"></span></label>
      </div>
      `;
    });
  }
  for (let s in g_sliders) {
    UpdateSlider(document.getElementById(s));
  }
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
  SetRuntimeConfig = Module.cwrap('SetRuntimeConfig', null, ['string']);
  SetExtraConfig = Module.cwrap('SetExtraConfig', null, ['string']);
  GetConfigString = Module.cwrap('GetConfigString', 'string', ['string']);
  GetConfigDouble = Module.cwrap('GetConfigDouble', 'number', ['string']);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = canvas.width / kScale;
  g_tmp_canvas.height = canvas.height / kScale;

  let keydown = function(e){
    if (e.key == ' ') {
      TogglePause();
      //e.preventDefault();
    }
    if(e.keyCode == 32 && e.target == document.body) {
      e.preventDefault();
    }
  };
  let keyup = function(e){
    if (e.key == ' ') {
      //e.preventDefault();
    }
  };
  let mouseclick = function(e){
    let x = e.offsetX / canvas.width;
    let y = 1 - e.offsetY / canvas.height;
    let r = GetConfigDouble("spawn_r");
    let extent = GetConfigDouble("extent");
    if (isNaN(r)) r = 0.05;
    if (isNaN(extent)) extent = 1;
    x *= extent;
    y *= extent;
    r *= extent;
    Spawn(x, y, r);
  };

  canvas.addEventListener('click', mouseclick, false);
  window.addEventListener('keydown', keydown, false);
  window.addEventListener('keyup', keyup, false);
  [
    input_conf,
  ].forEach(b => {
    b.addEventListener('keydown', function(e){e.stopPropagation();}, false);
    b.addEventListener('keyup', function(e){e.stopPropagation();}, false);
  });
  [
    window.button_pause, window.button_16, window.button_32,
    window.button_64, window.button_128, 
  ].forEach(b => {
    b.addEventListener('keydown', function(e){
      if (e.key == ' ') {
        e.preventDefault();
      }
    }, false);
    b.addEventListener('keyup', function(e){
      if (e.key == ' ') {
        e.preventDefault();
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
