#!/bin/sh

a=/tmp/a.$$
b=/tmp/b.$$
trap 'rm $a $b; exit 1' 1 2 3 4 15

ls /u/fall/*0.vtk > $a
ls /u/fall/*0.csv > $b

paste $a $b | awk '
{ 
  o = $1
  sub(/.*\//, "", o)
  cmd = sprintf("./rad %s %s > %s", q($1), q($2), q(o))
  r = system(cmd)
  print cmd
  if (r != 0) {
     printf "failed %s\n", cmd | "cat >&2"
     exit(1)
  }
}

function q(s)
{
        return "\"" s "\""
}
'

rm $a $b
