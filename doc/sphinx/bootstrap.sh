#!/bin/sh

awk '
BEGIN {
    printf "IMG_PDF = \\\n"
}

sub(/^\.\..* image:: /, "") {
    sub(/\.[^.]*$/, "")
    dirname = FILENAME
    sub(/\/[^\/]*$/, "", dirname)
    input = dirname "/" $0 ".svg"
    cmd = "test -r '\''" input "'\''"
    if (system(cmd) == 0) {
        output = dirname "/" $0 ".pdf"
        print output " \\" | "sort"
    }
}

END {
    close("sort")
    printf "\n"
}

' `find src -name '*.rst'` > img_pdf.mk

echo >&2 img_pdf.mk

IFS=,
set --- "$@"
