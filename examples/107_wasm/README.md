# WebAssembly Examples

Examples compiled to [WebAssembly](https://webassembly.org/)
and run directly in the web browser.

* [explorer](https://cselab.github.io/aphros/wasm/aphros.html):
Full solver with custom configuration

* [electrochem](https://cselab.github.io/aphros/wasm/electrochem.html):
Water electrolysis with dissolved gases nucleating into bubbles

* [hydro](https://cselab.github.io/aphros/wasm/hydro.html):
Interacting liquid drops with or without coalescence

* [diffusion](https://cselab.github.io/aphros/wasm/diffusion.html):
Advection-diffusion equation

## Build and run

Uses
[emscripten](https://emscripten.org/docs/getting_started/downloads.html#sdk-download-and-install),
commands `emcc` and `emar` should be available.

First, build and install WASM library in `src/`

```
cd src
make js
```

Then build and open examples in a web browser:

```
make
make run
```
