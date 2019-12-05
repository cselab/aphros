#!/usr/bin/awk -f

NR>1  {
  s[$1] +=$6 / $5
  #s[$1] ++
  n[$1] ++
  q[$1] += $a
  }

END {
  for (i in s)
    print i, s[i] / (n[i]+0.001) , n[i], q[i]
  }
