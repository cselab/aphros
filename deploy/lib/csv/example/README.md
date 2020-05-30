# Intro

   split -p prefix -f field [csv..]

# Example

   ./split -f cl -p a.%.csv /u/pasc_nn096/traj_010[0-9].csv
   ./csv2sph -p '%.vtk' data/traj_0100.csv
   ./csv2ellipsoid -p '%.vtk' data/traj_0100.csv
