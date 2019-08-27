module main
  use iso_c_binding
  implicit none
  real(8), parameter :: CONST_TINY = 1d-25
  real(8), parameter :: CONST_PI = 3.14159265358979323846d0
contains
  function cbrt(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: cbrt
    cbrt = x**(1.0/3.0)
  end function cbrt
  function calc_v0(alpha, vma, vmb, vmc) result(v)
    implicit none
    real(8), intent(in) :: alpha, vma, vmb, vmc
    real(8) :: v, a, vm1, vm2, vm3, vm12
    a = min(alpha, 1d0 - alpha)
    v = 0d0
    if (a > 0d0) then
       vm1 = min(vma, vmb, vmc)
       vm3 = max(vma, vmb, vmc)
       vm2 = abs(1d0 - vm3 - vm1)
       vm12 = vm1 + vm2
       if (a < vm1) then
          v = a ** 3 / (6d0 * vm1 * vm2 * vm3)
       else if (a < vm2) then
          v = a * (a - vm1) / (2d0 * vm2 * vm3) + &
               vm1 ** 2 / (6d0 * vm2 * vm3 + CONST_TINY)
       else if (a < min(vm12, vm3)) then
          v = (a ** 2 * (3d0 * vm12 - a) + vm1 ** 2 * (vm1 - 3d0 * a) + &
               vm2 ** 2 * (vm2 - 3d0 * a)) / (6d0 * vm1 * vm2 * vm3)
       else if (vm3 < vm12) then
          v = (a ** 2 * (3d0 - 2d0 * a) + vm1 ** 2 * (vm1 - 3d0 * a) + &
               vm2 ** 2 * (vm2 - 3d0 * a) + vm3 ** 2 * (vm3 - 3d0 * a)) / &
               (6d0 * vm1 * vm2 * vm3)
       else
          v = (a - 0.5d0 * vm12) / vm3
       end if
    end if
    if (alpha > 0.5d0) v = 1d0 - v
  end function calc_v0

  function calc_alpha0(v, vma, vmb, vmc) result(alpha)
    implicit none
    real(8), intent(in) :: v, vma, vmb, vmc
    real(8) :: alpha, w, vm1, vm2, vm3, vm12, v1, v3, a0, a1, a2, q0, sp, th
    w = min(v, 1d0 - v)
    vm1 = min(vma, vmb, vmc)
    vm3 = max(vma, vmb, vmc)
    vm2 = abs(1d0 - vm3 - vm1)
    vm12 = vm1 + vm2
    v1 = vm1 ** 2 / (6d0 * vm2 * vm3 + CONST_TINY)
    if (w < v1) then
       alpha = cbrt(6d0 * vm1 * vm2 * vm3 * w)
    else if (w < v1 + (vm2 - vm1) / (2.0 * vm3)) then
       alpha = 0.5d0 * (vm1 + sqrt(vm1 ** 2 + 8d0 * vm2 * vm3 * (w - v1)))
    else
       alpha = 0d0
       if (vm3 < vm12) then
          v3 = (vm3 ** 2 * (3d0 * vm12 - vm3) + vm1 ** 2 * (vm1 - 3d0 * vm3) + &
               vm2 ** 2 * (vm2 - 3d0 * vm3)) / (6d0 * vm1 * vm2 * vm3)
       else
          v3 = 0.5d0 * vm12 / vm3
          if (v3 <= w) alpha = vm3 * w + 0.5d0 * vm12
       end if
       if (alpha == 0.0) then
          if (w < v3) then
             a2 = -3d0 * vm12
             a1 = 3d0 * (vm1 ** 2 + vm2 ** 2)
             a0 = -(vm1 ** 3 + vm2 ** 3 - 6d0 * vm1 * vm2 * vm3 * w)
          else
             a2 = -1.5d0
             a1 = 1.5d0 * (vm1 ** 2 + vm2 ** 2 + vm3 ** 2)
             a0 = -0.5d0 * (vm1 ** 3 + vm2 ** 3 + vm3 ** 3 - 6d0 * vm1 * vm2 * vm3 * w)
          end if
          q0 = (1d0/6d0) * (a1 * a2 - 3d0 * a0) - a2 ** 3 * (1d0/27d0)
          sp = sqrt(-1d0/3d0 * a1 + 1d0/9d0 * a2 ** 2)
          th = 1d0/3d0 * acos(q0 / (sp ** 3))
          alpha = 2d0 * sp * cos(th + (4d0/3d0 * CONST_PI)) - (1d0/3d0) * a2
       end if
    end if
    if (v > 0.5d0) alpha = 1d0 - alpha
  end function calc_alpha0

  function calc_v(alpha, a, b, c) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: alpha, a, b, c
    real(c_double) :: calc_v
    calc_v = calc_v0(alpha, a, b, c)
  end function calc_v

  function calc_alpha(v, a, b, c) bind(c)
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: v, a, b, c
    real(c_double) :: calc_alpha
    calc_alpha = calc_alpha0(v, a, b, c)
  end function calc_alpha
end module main


