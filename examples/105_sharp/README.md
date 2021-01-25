```
cmake .
make
make install
make -C data/text
python3 -m pip install --user opencv-python
tools/sharppng data/text/blur.png o.png
```
