__kernel void iter2(
    int start, int lead_y, int lead_z, __global const Scal* u,
    __global Scal* output) {
  const unsigned int ix = get_global_id(0);
  const unsigned int iy = get_global_id(1);
  const unsigned int iz = get_global_id(2);
  const unsigned int i = start + iz * lead_z + iy * lead_y + ix;
  output[i] = u[i] + 1;
}
