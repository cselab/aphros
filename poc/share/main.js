function encodeSafe(base) {
  return base
          .replace(/\+/g, '-')
          .replace(/\//g, '_');
}
function decodeSafe(base) {
  return base
          .replace(/-/g, '+')
          .replace(/_/g, '/');
}
function update() {
  let url = new URL(location.href);
  let config = url.searchParams.get('config');
  let info = "";
  let compressed;
  let text_input = document.getElementById('text_input');
  if (config) {
    info += "Input from ?config=<br>";
    compressed = config;
  } else {
    info += "Input from text area<br>";
    compressed = encodeSafe(LZString.compressToBase64(text_input.value));
  }
  document.getElementById('text_compressed').value = compressed;
  sample = LZString.decompressFromBase64(decodeSafe(compressed));
  info += "Size of compressed: " + compressed.length + "<br>";
  info += "Size of decompressed: " + sample.length + "<br>";
  info += "Decompressed:<br><pre>" + sample + "</pre><br>";
  url.search = "?config=" + compressed;
  info += "<a href='" + url.toString() + "'>Link</a>";
  document.getElementById("div_info").innerHTML = info;
}
function request(url) {
  let div = document.getElementById("div_info");
  function write(msg) {
    div.innerHTML += msg + "<br>";
  };

  let xhr = new XMLHttpRequest();

  xhr.open('GET', url);

  xhr.send();

  xhr.onload = function() {
    if (xhr.status != 200) {
      write(`Error ${xhr.status}: ${xhr.statusText}`);
    } else {
      write(`Done, got ${xhr.response.length} bytes`);
      write(xhr.response);
    }
  };

  xhr.onerror = function() {
    write("Request failed");
  };
}
function requestGithub() {
  request('https://raw.githubusercontent.com/cselab/aphros/master/Dockerfile');
}
function requestTinyurl() {
  request('https://tinyurl.com/api-create.php?url=google.com');
}
function clearURL() {
  let url = new URL(location.href);
  url.search = "";
  history.replaceState(null, null, url.toString());
}
