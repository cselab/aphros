__kernel void iter2(
    int start, int lead_y, int lead_z, __global const Scal* fcp,
    __global const Scal* fcr, Scal dot_r, Scal dot_r_prev,
    __global Scal* output) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int iz = get_global_id(2);
  const unsigned int i = start + iz * lead_z + iy * lead_y + ix;
  output[i] = fcr[i] + (dot_r / (dot_r_prev + 1e-100)) * fcp[i];
}
