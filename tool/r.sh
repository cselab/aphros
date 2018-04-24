#!/bin/bash

set -eu

items='
s/Vars\.h/vars.h/gI
'

for i in $items ; do
  s="\"$i\""
  echo "$s"
  l=`cat l`
  for f in $l ; do
    cmd="sed -i $s $f"
    eval "$cmd"
  done
done
