#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <emscripten.h>
#include <emscripten/html5.h>
#include <cmath>
#include <iostream>

const int nx = 640;
const int ny = 480;

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

int g_time = 0;
float g_speed = 10;
float g_angle = 0;

static int sq(int a) {
  return a * a;
}

static float clamp(float f) {
  return f < 0 ? 0 : f > 1 ? 1 : f;
}

extern "C" {

float MulSpeed(float factor) {
  return g_speed *= factor;
}
}

static void main_loop() {
  memset(screen, 0, nx * ny * 4);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      int xx = x + sin(g_angle) * nx / 3 - nx / 2;
      int yy = y + cos(g_angle) * ny / 3 - ny / 2;
      int r = 50;
      int d = sq(xx) + sq(yy) - sq(r);
      int w = sq(r) / 5;
      float f = clamp((w - float(d)) / w);
      char cr = (1 + sin(g_time * 0.15)) * 127 * f;
      char cg = (1 + sin(g_time * 0.2)) * 127 * f;
      char cb = (1 + sin(g_time * 0.3)) * 127 * f;
      uint32_t c = 0xff000000 | cb | (cg << 8) | (cr << 16);
      screen[nx * y + x] = c;
    }
  }
  if (g_time % 10 == 0) {
    std::cout << "t=" << g_time << " speed=" << g_speed << std::endl;
  }
  g_angle += g_speed * 1e-3;
  ++g_time;
  Copy_ToCanvas(screen, nx, ny);
}

int main() {
  emscripten_set_canvas_element_size("#canvas", nx, ny);
  emscripten_set_main_loop(main_loop, 10, 1);
}
