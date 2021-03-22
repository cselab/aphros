A tool to make grayscale images sharper.

```
cmake . && make && make install # or make -f Makefile_legacy install
python3 -m pip install --user opencv-python
(cd data/text && make)
tools/sharp data/text/blur.pgm o.pgm --csv_out o.csv
```
