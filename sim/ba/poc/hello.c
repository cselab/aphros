#include <stdio.h>
#include <vof.h>

int
main(void)
{
    double f;
    f = vof_cylinder(0.11229068016121002, -0.7670366859609534, 0.22461479590800626, 0.5, 0.07159605407756882, -0.8070152387956918, 0.586174384796926);
    printf("f: %g\n", f);
}
