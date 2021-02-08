#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <emscripten.h>
#include <emscripten/html5.h>

#define nx 640
#define ny 480

uint32_t screen[nx * ny];

void Copy_ToCanvas(uint32_t* ptr, int w, int h) {
  EM_ASM_(
      {
        let data = Module.HEAPU8.slice($0, $0 + $1 * $2 * 4);
        let context = Module['canvas'].getContext('2d');
        let imageData = context.getImageData(0, 0, $1, $2);
        imageData.data.set(data);
        context.putImageData(imageData, 0, 0);
      },
      ptr, w, h);
}

int t = 0;

static int sq(int a) {
  return a * a;
}

static void main_loop() {
  memset(screen, 0, nx * ny * 4); // Clear screen
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      int xx = x + sin(t * 0.1) * 150 - nx / 2;
      int yy = y + cos(t * 0.1) * 100 - ny / 2;
      char cr = 127 + sin(t * 0.15) * 127;
      char cg = 127 + sin(t * 0.2) * 127;
      char cb = 127 + sin(t * 0.3) * 127;
      const int r = 50;
      screen[nx * y + x] = (sq(xx) + sq(yy) < sq(r) ? 0xff000000 + cb + (cg << 8) + (cr << 16) : 0);
    }
  }
  ++t;
  Copy_ToCanvas(screen, nx, ny);
}

int main() {
  EMSCRIPTEN_RESULT r = emscripten_set_canvas_element_size("#canvas", nx, ny);
  emscripten_set_main_loop(main_loop, 30, 1);
  assert(r == EMSCRIPTEN_RESULT_SUCCESS);
  return 0;
}
