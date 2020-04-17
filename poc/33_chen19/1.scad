n = 1;
s = 200; // spacing
h = 210;
d1 = 51.3;
d2 = 129.0;

$fn = 10;
minkowski() {
difference() {
	     scale([s * n, s * n, h]) cube(1);
	     for (i=[0 : n - 1])
		 for (j=[0 : n - 1])
		     translate([s * (i + 0.5), s * (j + 0.5)]) cylinder(r1 = d1/2, r2 = d2/2, h = h);
	}
sphere(5);
}
