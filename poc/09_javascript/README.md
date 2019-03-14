# Intro

`*.in.js` and `*.in.html` files are preprocessed with `m4`. They
should include line `include(mh.m4)dnl`.


# Installl

	make install
	make test

# Plan

    function Parstr(nh, hp, eta) {
		n, nh, hp, eta, a0, t0, segmetns, [storage for particles]

		this.start = function(ends, i, j, a, t) { }
		this.step = function(n) { }
		this.converge = function(nmax, eps) { }
		this.eps = function() { }
		this.a = function() { }
		this.t = function() { }
		this.force = function () { }
		this.particle = function () { }
    }
