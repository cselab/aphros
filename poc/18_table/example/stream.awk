#!/usr/bin/awk -f

/^[\t ]*put/ {
    x = $2
    y = $3
    t[x] = y
}

/^[\t ]*get/ {
    x = $2
    print t[x]
}

/^[\t ]*remove/ {
    x = $2
    delete t[x]
}

/^[\t ]*length/ {
    n = 0
    for (i in t)
        n++
    print n
}
