    cl=1367 ; ./traj_recover.py $cl ./traj_*.csv > trajr_$cl.csv

    vtkin=../../sm_%04d.vtk ./vtkcutcl ../trajr_7850.csv
