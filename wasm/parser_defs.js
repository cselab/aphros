var SetConfig;
var GetConfig;
var input_conf = document.getElementById('input_conf');
var output = document.getElementById('output');
var outputerr = document.getElementById('outputerr');
var textarea_updated = true;
var text_shorturl = document.getElementById('text_shorturl');
var g_fullurl;
var last_set_config = 0;

function Draw() {
  if (textarea_updated) {
    UpdateFullUrl();
    res = SetConfig(input_conf.value);
    textarea_updated = false;
    if (output) {
      output.value = GetConfig();
      if (res == 0 && last_set_config) {
        PrintError("OK");
      }
    }
    last_set_config = res;
  }
}

function ClearOutput() {
  output.value = '';
  outputerr.value = '';
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
  history.replaceState(null, null, fullurl);
  if (g_fullurl != fullurl) {
    ClearShortUrl();
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
function ClearShortUrl() {
  text_shorturl.value = "";
}
function UpdateShortUrl() {
  let fullurl = GetFullUrl(input_conf.value);
  g_fullurl = fullurl;
  GetShortUrl(fullurl, function(response) {
    text_shorturl.value = response;
  });
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
function Compress(text) {
  return EncodeSafe(LZString.compressToBase64(text));
}

function Decompress(compressed) {
  return LZString.decompressFromBase64(DecodeSafe(compressed));
}
function PostRun() {
  SetConfig = Module.cwrap('SetConfig', 'int', ['string']);
  GetConfig = Module.cwrap('GetConfig', 'string', []);

  let url = new URL(location.href);
  let compressed = url.searchParams.get('config');
  if (compressed) {
    input_conf.value = Decompress(compressed);
  }
  UpdateFullUrl();
}

var Module = {
  preRun: [],
  postRun: [PostRun],
  printErr: (function(text) {
    ClearOutput();
    return function(text) {
      if (arguments.length > 1) {
        text = Array.prototype.slice.call(arguments).join(' ');
      }
      PrintError(text);
    };
  })(),
  setStatus: function(text) {},
};
