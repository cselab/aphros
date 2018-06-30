#!/bin/bash
alias gerris2D=ge.run2d

if test x$donotrun != xtrue; then
    for La in 120 1200 12000; do
	tmax=`echo $La | awk '{print 0.8*0.8/sqrt(0.8/$1)}'`
	if sed "s/end = TMAX/iend = 1/g" < $1 | gerris2D -DLEVEL=5 -DLAPLACE=$La -DDT=0 - |\
           sed "s/iend = 1/end = $tmax/" | gerris2D - > /dev/null; then :
	else
	    exit 1
	fi
    done

    La=12000
    for level in 3 4 6 7; do
	tmax=`echo $La | awk '{print 0.8*0.8/sqrt(0.8/$1)}'`
	if sed "s/end = TMAX/iend = 1/g" < $1 | gerris2D -DLEVEL=$level -DLAPLACE=$La -DDT=1e-9 - |\
           sed "s/iend = 1/end = $tmax/" | gerris2D - > sim-$level; then : 
	else
	    exit 1
	fi
    done
fi

rm -f convergence kconvergence
La=12000
for level in 3 4 5 6 7; do
    if awk -v level=$level < E-$La-$level '{
             max2 = $3
             maxi = $4
           }END{print 0.8*2**level, max2, maxi}' >> convergence; then : 
    else
	exit 1
    fi
    if awk -v level=$level < EK-$La-$level '{
             max2 = $3
             maxi = $4
           }END{print 0.8*2**level, max2/2.5, maxi/2.5}' >> kconvergence; then : 
    else
	exit 1
    fi
done

if cat <<EOF | gnuplot ; then :
    set term postscript eps color lw 3 solid 20
    set output 'laplace.eps'
    set xlabel 'Tau'
    set ylabel 'U(D/sigma)^1/2'
    set logscale y
    plot 'La-120-5' w l t "La=120", 'La-1200-5' w l t "La=1200", 'La-12000-5' w l t "La=12000"
    set output 'curvature.eps'
    set ylabel 'Curvature standard deviation'
    plot [][1e-12:]'K-120-5' u 1:4 w l t "La=120", 'K-1200-5' u 1:4 w l t "La=1200", 'K-12000-5' u 1:4 w l t "La=12000"
    set output 'convergence.eps'
    set xlabel 'D'
    set ylabel 'Shape error'
    set logscale x
    set xtics 2
    plot [5:120]'convergence' u 1:2 w lp t "RMS" ps 3, 'convergence' u 1:3 w lp t "Max" ps 3, 0.2/(x*x) t "Second order"
    set output 'kconvergence.eps'
    set ylabel 'Relative curvature error'
    set logscale x
    plot [5:120]'kconvergence' u 1:3 w lp t "" ps 3, 0.6/(x*x) t "Second order"
EOF
else
    exit 1
fi

for f in La-120-5 La-1200-5 La-12000-5; do
    if tail -n 1 < $f | awk -v tolerance=$2 '{if ($2 > tolerance) exit (1);}'; then :
    else
	exit 1
    fi
done

if cat <<EOF | python ; then :
from check import *
from sys import *
if (Curve('convergence',1,3) - Curve('convergence.ref',1,3)).max() > 1e-6:
    exit(1)
if (Curve('kconvergence',1,3) - Curve('kconvergence.ref',1,3)).max() > 1e-6:
    exit(1)
EOF
else
   exit 1
fi
