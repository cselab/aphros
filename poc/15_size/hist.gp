plot [0.01:][] \
     "data/slip.gyr" w lp, \
     "data/fall.gyr" w lp, \
     0.02*(x/0.02)**(-7/3.0), \
     0.01*(x/0.02)**(-10/3.0)
