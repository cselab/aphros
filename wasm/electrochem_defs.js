var g_lines;
var g_lines_ptr;
var SetRuntimeConfig;
var GetLines;
var Spawn;
var g_tmp_canvas;
var kScale = 1;
var kMarginX = 32; // XXX must match g_canvas->xmargin

function GetExtraConfig() {
  return `
set double rho1 1
set double rho2 0.1

set int nsteps 2

set double cfl 0.7

set double growth_rate 40

set double resist1 1
set double resist2 100

set int visinterp 1

set double surface_height 0.8

set double visvort 0.01
set double vispot 0
set double vistracer 2
set double visvf 0.7
set double visnucl 1

set double cmax 0.6
set double nucl_noise 0.01
set double nucl_dt 0.01
set double nucl_c_add 0.5
set double nucl_c_remove 0.3
set double nucl_vf_remove 0.8


set double hypre_symm_tol 0.1
set int hypre_symm_maxiter 50
set int hypre_symm_miniter 10
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
  cc.push(SetSigma(window.range_sigma.value));
  cc.push(SetMu(window.range_mu.value));
  cc.push(SetRate(window.range_rate.value));
  cc.push(SetDiffusion(window.range_diffusion.value));
  cc.push(SetGravity(window.range_gravity.value));
  ResetButtons();
  let button = document.getElementById('button_' + nx);
  if (button) {
    button.className = "button pressed";
  }
  conf += cc.join('');
  SetExtraConfig(conf);
  Module.ccall('SetMesh', null, ['number'], [nx])
}

function DrawBackground() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  { // electrodes
    let w = 30;
    let lw = 2;
    let lwh = lw / 2;
    class Text {
      constructor(text, size, dx, dy) {
        this.text = text;
        this.size = size;
        this.dx = dx;
        this.dy = dy;
      }
    };
    let electrode = function(rect, texts) {
      ctx.lineWidth = lw;
      ctx.strokeStyle="black";
      ctx.fillStyle="#FFC61E";

      ctx.fillRect(...rect);
      ctx.strokeRect(...rect);

      ctx.fillStyle = 'black';
      for (let i = 0; i < texts.length; ++i) {
        let t = texts[i];
        ctx.font = t.size + 'px  Arial, Helvetica, sans-serif';
        ctx.fillText(t.text, rect[0] + t.dx, rect[1] + rect[3] / 2 + t.dy);
      }
    };
    electrode([kMarginX - w - lwh, lwh, w, canvas.height - lw],
        [
          new Text("\u2013", 30, 6, 0), // minus
          new Text("H\u2082", 22, 2, 50), // hydrogen
        ]);
    electrode([canvas.width - kMarginX + lwh, lwh, w, canvas.height - lw],
        [
          new Text("+", 30, 6, 0), // plus
          new Text("O\u2082", 22, 2, 50), // oxygen
        ]);
  }
}

function Draw() {
  let canvas = Module['canvas'];
  let ctx = canvas.getContext('2d');
  ctx.drawImage(g_tmp_canvas, kMarginX, 0,
      canvas.width - 2 * kMarginX, canvas.height);

  {
    let lw = 2;
    let lwh = lw / 2;
    ctx.lineWidth = lw;
    ctx.strokeStyle="black";
    ctx.strokeRect(kMarginX - lwh, lwh,
        canvas.width - 2 * kMarginX + lw, canvas.height - lw);
  }

  g_lines = new Uint16Array(Module.HEAPU8.buffer, g_lines_ptr, g_lines_max_size);
  let size = GetLines(g_lines.byteOffset, g_lines.length);
  ctx.lineWidth = 2;
  ctx.strokeStyle = "black";
  for (let i = 0; i + 3 < size; i += 4) {
    ctx.beginPath();
    ctx.moveTo(g_lines[i], g_lines[i + 1]);
    ctx.lineTo(g_lines[i + 2], g_lines[i + 3]);
    ctx.stroke();
  }
}

function PostRun() {
  g_lines_max_size = 1024 * 100;
  g_lines_ptr = Module._malloc(g_lines_max_size * 2);
  Spawn = Module.cwrap('Spawn', null, ['number', 'number', 'number']);
  GetLines = Module.cwrap('GetLines', 'number', ['number', 'number']);

  let canvas = Module['canvas'];
  g_tmp_canvas = document.createElement('canvas');
  g_tmp_canvas.width = (canvas.width - 2 * kMarginX) / kScale;
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
    let x = (ev.offsetX - kMarginX) / (canvas.width - 2 * kMarginX);
    let y = 1 - ev.offsetY / canvas.height;
    Spawn(x, y, 0.05);
  };

  window.addEventListener('keydown', keydown, false);
  window.addEventListener('keyup', keyup, false);
  canvas.addEventListener('click', mouseclick, false);
  Init(32);
  DrawBackground();
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
