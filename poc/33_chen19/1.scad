n = 5;
d = 300;
h = 210;

r1 = 51.3;
r2 = 129.0;

difference() {
             scale([d * n, d * n, h]) cube(1);
	     for (i=[0 : n - 1])
	      	 for (j=[0 : n - 1])
	     	     translate([d * (i + 0.5), d * (j + 0.5), 0]) cylinder(r1 = r1, r2 = r2, h = h);
	}
