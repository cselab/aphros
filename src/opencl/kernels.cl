__kernel void inner_to_buf_dim3(
    int start, int lead_y, int lead_z, __global const Scal* field,
    __global Scal* buf) {
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int nx = get_global_size(0);
  const int ny = get_global_size(1);
  const int nz = get_global_size(2);
  int i = 0;

#define X(IX, IY, IZ, NX, NY, NZ)                                          \
  if (IX == 0) {                                                           \
    buf[i + NY * IZ + IY] = field[start + iz * lead_z + iy * lead_y + ix]; \
  }                                                                        \
  i += NY * NZ;                                                            \
  if (IX == NX - 1) {                                                      \
    buf[i + NY * IZ + IY] = field[start + iz * lead_z + iy * lead_y + ix]; \
  }                                                                        \
  i += NY * NZ;

  X(ix, iy, iz, nx, ny, nz);
  X(iy, iz, ix, ny, nz, nx);
#if DIM > 2
  X(iz, ix, iy, nz, nx, ny);
#endif
#undef X
}

__kernel void buf_to_halo_dim3(
    int start, int lead_y, int lead_z, __global const Scal* buf,
    __global Scal* field) {
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);
  const int nx = get_global_size(0);
  const int ny = get_global_size(1);
  const int nz = get_global_size(2);
  int i = 0;

#define X(IX, IY, IZ, NX, NY, NZ)                                          \
  if (IX == 0) {                                                           \
    const int t = IX - 1;                                                  \
    const int IX = t;                                                      \
    field[start + iz * lead_z + iy * lead_y + ix] = buf[i + NY * IZ + IY]; \
  }                                                                        \
  i += NY * NZ;                                                            \
  if (IX == NX - 1) {                                                      \
    const int t = IX + 1;                                                  \
    const int IX = t;                                                      \
    field[start + iz * lead_z + iy * lead_y + ix] = buf[i + NY * IZ + IY]; \
  }                                                                        \
  i += NY * NZ;

  X(ix, iy, iz, nx, ny, nz);
  X(iy, iz, ix, ny, nz, nx);
#if DIM > 2
  X(iz, ix, iy, nz, nx, ny);
#endif
#undef X
}

__kernel void field_max(
    int start, int lead_y, int lead_z, __global const Scal* u,
    __global Scal* output) {
  const int ngx = get_num_groups(0);
  const int ngy = get_num_groups(1);
  const int gx = get_group_id(0);
  const int gy = get_group_id(1);
  const int gz = get_group_id(2);
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int i = start + iz * lead_z + iy * lead_y + ix;
  output[gz * ngy * ngx + gy * ngx + gx] = work_group_reduce_max(u[i]);
}

__kernel void field_sum(
    int start, int lead_y, int lead_z, __global const Scal* u,
    __global Scal* output) {
  const int ngx = get_num_groups(0);
  const int ngy = get_num_groups(1);
  const int gx = get_group_id(0);
  const int gy = get_group_id(1);
  const int gz = get_group_id(2);
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int i = start + iz * lead_z + iy * lead_y + ix;
  output[gz * ngy * ngx + gy * ngx + gx] = work_group_reduce_add(u[i]);
}

__kernel void field_dot(
    int start, int lead_y, int lead_z, __global const Scal* u,
    __global const Scal* v, __global Scal* output) {
  const int ngx = get_num_groups(0);
  const int ngy = get_num_groups(1);
  const int gx = get_group_id(0);
  const int gy = get_group_id(1);
  const int gz = get_group_id(2);
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int i = start + iz * lead_z + iy * lead_y + ix;
  output[gz * ngy * ngx + gy * ngx + gx] = work_group_reduce_add(u[i] * v[i]);
}
