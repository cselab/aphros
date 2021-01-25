The Voronoi Implicit Interface Method for computing multiphase physics
Robert I. Sayea and James A. Sethiana
2011

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3241776/

https://www.pnas.org/content/pnas/108/49/19498/F6.large.jpg


# Workflow

* download image from pnas 

      ./0_wget

* crop image to create `a0.png` etc

      ./1_crop

* open `voronoi/online/voronoi.html` in browser

* make the browser window transparent and put `a0.png` below

* create points with mouse click and move them with shift+click,
  adjust to match the image

* copy the sites from the edit box in __Input data__ to file `voronoi/online/sites0`

* parse the file to extract coordinates of sites (from `voronoi/online/`)

      ./sites2points sites0 > points0

* create Voronoi diagram of points (from `voronoi/`)

    ./voronoi.py online/points0

* run simulation in `sim` reading field from `ref/voronoi/cl.dat`

* read reference colors (from `ref`)

      ./points2colors voronoi/online/points0 a0.png > colors0
      ./points2colors voronoi/online/points1 a0.png 3 > colors1

* plot simulation output with reference colors (from `sim/vispy/`)

      ./plot.py sm_0200.vtk ../../ref/colors0


