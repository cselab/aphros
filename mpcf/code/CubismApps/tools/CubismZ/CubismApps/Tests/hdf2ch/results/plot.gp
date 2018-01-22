#data=$1
#x=$2
#y=$3

print data

#set terminal jpeg 
#set output data.".jpg"

set terminal pdf
set output data.".pdf"

set yrange [0:240]
set ytics 0,20,240
#set xrange [1:10000]
#set xtics 1,10,10000
set xrange [1:1000]
set xtics 1,10,1000
set grid xtics ytics

set logscale x
#set logscale y
set xlabel xlab
set ylabel ylab
set title head

psize=0.2
lsize=5

plot data."_wavz_res.txt" u x:y w lp ps psize lw lsize title 'WAVZ-3', \
     data."_zfp_res.txt" u x:y w lp ps psize lw lsize title 'ZFP-0.5.0', \
     data."_sz1_res.txt" u x:y w lp ps psize lw lsize title 'SZ-1.4.8b', \
     data."_fpzip1_res.txt" u x:y w lp ps psize lw lsize title 'FPZIP-1.1.0', \
     data."_wavz1_res.txt" u x:y w lp ps psize lw lsize title 'WAVZ-1', \

#     data."_sz_res.txt" u x:y w lp ps psize lw lsize title 'SZ-0.5.14', \
#     data."_fpzip_res.txt" u x:y w lp title 'FPZIP-1.0.1', \



#  plot data."_wavz_res.txt" u x:y w lp title 'WAVZ'
#replot data."_fpzip_res.txt" u x:y w lp title 'FPZIP'
#replot data."_zfp_res.txt" u x:y w lp title 'ZFP'
#replot data."_sz_res.txt" u x:y w lp title 'SZ'
