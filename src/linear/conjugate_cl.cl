__kernel void iter2(
    int start, int lead_y, int lead_z, __global const Scal* fcp,
    __global const Scal* fcr, Scal dot_r, Scal dot_r_prev,
    __global Scal* output) {
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int i = start + iz * lead_z + iy * lead_y + ix;
  output[i] = fcr[i] + (dot_r / (dot_r_prev + 1e-100)) * fcp[i];
}

__kernel void linear(
    int start, int lead_y, int lead_z, __global const Scal* fc_system,
    __global const Scal* fcu, __global Scal* output, Scal k0, Scal k1) {
  const int ix = get_global_id(0);
  const int iy = get_global_id(1);
  const int iz = get_global_id(2);
  const int i = start + iz * lead_z + iy * lead_y + ix;
  const int ixm = i - 1;
  const int ixp = i + 1;
  const int iym = i - lead_y;
  const int iyp = i + lead_y;
  const int izm = i - lead_z;
  const int izp = i + lead_z;
  const int k = i * (DIM * 2 + 2);
  output[i] = k0 * (fc_system[k] * fcu[i] + //
                    fc_system[k + 1] * fcu[ixm] + fc_system[k + 2] * fcu[ixp] +
                    fc_system[k + 3] * fcu[iym] + fc_system[k + 4] * fcu[iyp] +
                    fc_system[k + 5] * fcu[izm] + fc_system[k + 6] * fcu[izp]) +
              k1 * fc_system[k + 7];
}
