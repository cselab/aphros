# split stat
  2 1
  32574 2
    985 3
     42 4
      5 5

# Commands

    ch.split -p a_%.csv -f cl /u/pasc_nn096/traj_*.csv
    cp a_0001.csv a_0000.csv
    ch.color -p a_%.vtk -k cl -f split a_*0.csv -- /u/pasc_nn096/sm_*0.vtk
    ./split a_0*.csv  > stat


    ./split ~/tmp/qh/a_0{0..2}??.csv > stathvof
    ./split ~/tmp/qhvof/a_0{0..2}??.csv > stathvof
    ./split ~/tmp/qhvof/a_0{0..2}??.csv > stathvof

# awk

    cawk 'prog' [file ..]