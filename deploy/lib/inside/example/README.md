Bonds
```
./plane -a 0 -n 0 1 0 ~/g.off | ./2line
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