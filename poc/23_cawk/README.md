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


    ./split ~/tmp/qh/a_0{0..2}??.csv > stath
    ./split ~/tmp/qhvof/a_0{0..2}??.csv > stathvof
    ./split ~/tmp/qhvof/a_0{0..2}??.csv > stathvof

    ch.filter -f cl -p sf_%.vtk 3598 12874 -- sm_0{158..163}.vtk

    co.gyration -p a_%.csv /u/pasc_nn096/traj_*.csv
    ch.color -p a_%.vtk -k cl -f r rg asphericity -- a_*0.csv -- /u/pasc_nn096/sm_*0.vtk

    ch.cawk '$r < 0.2 {print $r, $vf/(ARGC - 1)}' /u/pasc_nn192/traj_0[2]*.csv  | ch.hist --range 0 0.2 --bins 50 > b

# awk

    cawk 'prog' [file ..]
