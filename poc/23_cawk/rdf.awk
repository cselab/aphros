#/bin/sh

me=rdf
dmin=0.01
dmax=0.1

for i
do
echo >&2 "$i"
ch.cawk -v me="$me" -v dmin=$dmin -v dmax=$dmax '
$r < 0.1 && $r > 1/192 {
    n += 0
    x[n] = $x
    y[n] = $y
    z[n] = $z
    r[n] = $r
    n++
}

END {
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            if (i == j)
                continue
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            d = dx*dx + dy*dy + dz*dz
            d = sqrt(d)
            if (d != 0 && d > dmin && d < dmax)
                print d, 1/d^2
        }
}

' "$i"
done
