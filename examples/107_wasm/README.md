# WebAssembly Examples

Examples compiled to [WebAssembly](https://webassembly.org/)
that directly in the web browser.

* `electrochem`: Water electrolysis with dissolved gases nucleating into bubbles
[Online demo](https://cselab.github.io/aphros/wasm/electrochem.html)

* `hydro`: Interacting liquid drops with or without coalescence
[Online demo](https://cselab.github.io/aphros/wasm/hydro.html)

* `diffusion`: Advection-diffusion equation
[Online demo](https://cselab.github.io/aphros/wasm/diffusion.html)


## Build and run

Uses
[emscripten](https://emscripten.org/docs/getting_started/downloads.html#sdk-download-and-install),
commands `emcc` and `emar` should be available. Build and open in a
web browser:

```
make
emrun --serve_after_exit hydro.html
```
