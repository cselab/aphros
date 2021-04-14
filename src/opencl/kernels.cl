__kernel void diffusion(
    int start, int lead_x, Scal k, //
    __global const Scal* u, __global Scal* output) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int i = start + iy * lead_x + ix;

  const Scal uc = u[i];
  const Scal uxm = u[start + iy * lead_x + (ix - 1)];
  const Scal uxp = u[start + iy * lead_x + (ix + 1)];
  const Scal uym = u[start + (iy - 1) * lead_x + ix];
  const Scal uyp = u[start + (iy + 1) * lead_x + ix];

  output[i] = uc + (uxm + uxp + uym + uyp - 4 * uc) * k;
}

__kernel void inner_to_buf(
    int start, int lead_x, __global const Scal* field, __global Scal* buf) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int nx = get_global_size(0);
  const unsigned int ny = get_global_size(1);
  int i = 0;
  if (iy == 0) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;
  if (iy == 1) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;
  if (iy == ny - 2) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;
  if (iy == ny - 1) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;

  if (ix == 0) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
  if (ix == 1) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
  if (ix == nx - 2) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
  if (ix == nx - 1) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
}

__kernel void inner_to_buf1(
    int start, int lead_x, __global const Scal* field, __global Scal* buf) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int nx = get_global_size(0);
  const unsigned int ny = get_global_size(1);
  int i = 0;
  if (iy == 0) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;
  if (iy == ny - 1) {
    buf[i + ix] = field[start + iy * lead_x + ix];
  }
  i += nx;

  if (ix == 0) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
  if (ix == nx - 1) {
    buf[i + iy] = field[start + iy * lead_x + ix];
  }
  i += ny;
}

__kernel void buf_to_halo1(
    int start, int lead_x, __global const Scal* buf, __global Scal* field) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int nx = get_global_size(0);
  const unsigned int ny = get_global_size(1);
  int i = 0;
  if (iy == 0) {
    field[start + (iy - 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;
  if (iy == ny - 2) {
    field[start + (iy + 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;

  if (ix == 0) {
    field[start + iy * lead_x + ix - 2] = buf[i + iy];
  }
  i += ny;
  if (ix == nx - 2) {
    field[start + iy * lead_x + ix + 2] = buf[i + iy];
  }
  i += ny;
}


__kernel void buf_to_halo(
    int start, int lead_x, __global const Scal* buf, __global Scal* field) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int nx = get_global_size(0);
  const unsigned int ny = get_global_size(1);
  int i = 0;
  if (iy == 0) {
    field[start + (iy - 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;
  if (iy == 1) {
    field[start + (iy - 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;
  if (iy == ny - 2) {
    field[start + (iy + 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;
  if (iy == ny - 1) {
    field[start + (iy + 2) * lead_x + ix] = buf[i + ix];
  }
  i += nx;

  if (ix == 0) {
    field[start + iy * lead_x + ix - 2] = buf[i + iy];
  }
  i += ny;
  if (ix == 1) {
    field[start + iy * lead_x + ix - 2] = buf[i + iy];
  }
  i += ny;
  if (ix == nx - 2) {
    field[start + iy * lead_x + ix + 2] = buf[i + iy];
  }
  i += ny;
  if (ix == nx - 1) {
    field[start + iy * lead_x + ix + 2] = buf[i + iy];
  }
  i += ny;
}

__kernel void field_sum(
    int start, int lead_x, __global const Scal* u, __global Scal* output) {
  const unsigned int ngx = get_num_groups(0);
  const unsigned int gx = get_group_id(0);
  const unsigned int gy = get_group_id(1);
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int i = start + iy * lead_x + ix;
  output[gy * ngx + gx] = work_group_reduce_add(u[i]);
}

__kernel void field_dot(
    int start, int lead_x, __global const Scal* u, __global const Scal* v,
    __global Scal* output) {
  const unsigned int ngx = get_num_groups(0);
  const unsigned int gx = get_group_id(0);
  const unsigned int gy = get_group_id(1);
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int i = start + iy * lead_x + ix;

  output[gy * ngx + gx] = work_group_reduce_add(u[i] * v[i]);
}
