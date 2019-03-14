# Intro

`*.in.js` and `*.in.html` files are preprocessed with `m4`. They
should include line `include(mh.m4)dnl`.


# Installl

	make install
	make test

# Plan

	function Parstr(nh, hp, eta) {
		segmetns, [storage for particles], n, nh, hp, eta
		a0, t0

		this.new  = function(a, t, ends, i, j) { }
		this.step = function(n) { }
		this.converge = function(nmax, eps) { }

		this.eps = function() { }
		this.a = function() { }
		this.t = function() { }
		this.force = function () { }
		this.particle = function () { }
	}
