n = 32;
s = 200; // spacing
h = 210;
d1 = 51.3;
d2 = 129.0;

linear_extrude(height = h)
difference() {
	     scale([s * n, s * n]) square(1);
	     for (i=[0 : n - 1])
		 for (j=[0 : n - 1])
		     translate([s * (i + 0.5), s * (j + 0.5)]) circle(r = d1/2);
	}
