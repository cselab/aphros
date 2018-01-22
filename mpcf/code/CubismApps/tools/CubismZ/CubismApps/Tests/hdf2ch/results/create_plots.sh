#!/bin/bash -l

gnuplot --persist -e "data='data_010000-p.h5'"      -e "head='p, 10k steps'"      -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-rho.h5'"    -e "head='rho, 10k steps'"    -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-a2.h5'"     -e "head='a2, 10k steps'"     -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-divU.h5'"   -e "head='divU, 10k steps'"   -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-E.h5'"      -e "head='E, 10k steps'"      -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-Omegax.h5'" -e "head='Omegax, 10k steps'" -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_010000-Ux.h5'"     -e "head='Ux, 10k steps'"     -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 


gnuplot --persist -e "data='data_005000-p.h5'"      -e "head='p, 5k steps'"      -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-rho.h5'"    -e "head='rho, 5k steps'"    -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-a2.h5'"     -e "head='a2, 5k steps'"     -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-divU.h5'"   -e "head='divU, 5k steps'"   -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-E.h5'"      -e "head='E, 5k steps'"      -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-Omegax.h5'" -e "head='Omegax, 5k steps'" -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
gnuplot --persist -e "data='data_005000-Ux.h5'"     -e "head='Ux, 5k steps'"     -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 

#gnuplot --persist -e "data='data-301-p.h5'" -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 
#gnuplot --persist -e "data='data-301-g.h5'" -e "x=2" -e "y=9" -e "xlab='Compression Ratio'" -e "ylab='PSNR'" plot.gp 

