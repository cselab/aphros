# WebAssembly Examples

Examples compiled to [WebAssembly](https://webassembly.org/)
and run directly in the web browser.

* [explorer](https://cselab.github.io/aphros/wasm/aphros.html):
Aphros Explorer: full solver with custom configuration.

* [electrochem](https://cselab.github.io/aphros/wasm/electrochem.html):
Water electrolysis with dissolved gases nucleating into bubbles.

* [hydro](https://cselab.github.io/aphros/wasm/hydro.html):
Interacting liquid drops with or without coalescence.

* [diffusion](https://cselab.github.io/aphros/wasm/diffusion.html):
Advection-diffusion equation.

* [parser](https://cselab.github.io/aphros/wasm/parser.html?config=M4UwLgBAlgdpUQAwChSQCYHsCuAjANiBOhAIwB0ArKuBAG4gDGkdZEATBAMw2TBgAnWAHMIwCGBAAPMLzGCRYgPoBbCACJkECPlghSWnXs4ASYMk1p5QmKOB4IJhCdZm5_G6JB0AhvkcAFAEmJABUHACUEXJYeIQQoihWDMwJAj50UGAAno6JQA=):
Parser of configuration files.

## Build and run

These examples use [emscripten](https://emscripten.org/docs/getting_started/downloads.html#sdk-download-and-install).

First, build and install `aphrosjs` library in `src/`

```
cd src
make js
```

Then build and open examples in a web browser:

```
make
make run
```
