Bonds
```
./plane -a 0 -n 0 1 0 ~/g.off | ./2line -o gnuplot > q
./plane -a 0 -n 0 1 0 ~/g.off | ./2off > q.off
```

Cones
```
./plane -a 2.5 -n 0 0 1 ~/g.off | ./2off > c.off
```

Big side tubes
```
./plane -a 0 -n 1 0 0 ~/g.off | ./2off > c.off
```

Generate XZ cross-section of the bonds
```
./gen data/01.line data/01.idx
```

Generate big side tubes
y = 0, z = 4, R = 2.5, lo = -7.03178e+00, hi = 2.334422e+01
