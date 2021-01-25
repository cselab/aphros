```
cmake .
make
make install
make -C data/text
python3 -m pip install --user h5py opencv-python
tools/sharppng data/text/blur.png o.png
```

# UZH

```
convert i.png -brightness-contrast -5x5 ii.png
cfl=0.5 steps=5 tools/sharppng ii.png o.png --vtk_out_march sm.vtk
```
