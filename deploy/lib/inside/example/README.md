Bonds
```
./plane -a 0 -n 0 1 0 ~/g.off | ./2line -o gnuplot > q
./plane -a 0 -n 0 1 0 ~/g.off | ./2off > q.off
gnuplot> plot "<./plane -a 0 -n 0 1 0 ~/g.off | ./2line -o gnuplot" u 1:3 w l
```

Cones
```
./plane -a 2.5 -n 0 0 1 ~/g.off | ./2off > c.off
```

Big side tubes
```
./plane -a 0 -n 1 0 0 ~/g.off | ./2off > c.off
```

Brick
x: -2.3344220000000000e+01 2.0758048711630039e+01
y: -5 5
z: -6.5 6.5

Generate XZ cross-section of the bonds
```
./gen data/01.line data/01.idx
```

Generate big side tubes
```
./prim2off < data/big.dat > q.off
```
