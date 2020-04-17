s = 200; // spacing
h = 210;
d1 = 51.3;
d2 = 129.0;

$fn = 30;
minkowski() {
	  translate([s * 0.5, s * 0.5]) cylinder(r1 = d1/2, r2 = d2/2, h = h);
	  sphere(10);
}
