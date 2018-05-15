plot \
"<paste ref/{t,dk}" w l lw 3 lc "b" t "ref" \
, "<paste sc64i2cfl025/{t,dk}" w l lw 2 t "64bi2025" \
, "<paste sc64q/{t,dk}" w l lw 2 t "64q" \
, "<paste sc/{t,dk}" w l lw 2 t "64c" \

#, "<paste sc64/{t,dk}" w l lw 2 t "64" \
#, "<paste sc64i2/{t,dk}" w l lw 2 t "64bi2" \
#, "<paste sc64i5/{t,dk}" w l lw 2 t "64bi5" \
#, "<paste sc/{t,dk}" w l lw 2 t "64b" \
#, "<paste ../sim07_256/sc/{t,dk}" w l lw 2 t "256" \
#, "<paste ../sim07_b1/sc/{t,dk}" w l lw 2 t "128,b1" \
#, "<paste sc64/{t,dk}" w l lw 3 t "64" \
#, "<paste ../sim07_512/sc256/{t,dk}" w l lw 2 t "256" \

#set xrange [0:2]
set yrange [-0.002:0.014]
