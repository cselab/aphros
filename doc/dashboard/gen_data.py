#!/usr/bin/env python3

import glob
import re
import subprocess
import collections
import argparse


def gen(wildcard, outpath, urlbase):
    dd = sorted(glob.glob(wildcard))
    print("Found {:} files by '{}'".format(len(dd), wildcard))

    res = []

    for d in dd:
        with open(d) as f:
            text = f.read()
        rev = re.findall("(?:hydro|aphros) (........)", text)[0]
        runtime = float(re.findall("all \[([\d.]*) s", text)[0])
        url = urlbase + rev
        p = subprocess.Popen(
            ["git", "log", "-1", "--pretty=%h %as", rev, "--"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out = str(p.stdout.read().decode("utf-8")).strip().split()
        if len(out):
            gitrev, gitdate = out
            if rev == gitrev:
                res.append({
                    'rev': rev,
                    'runtime': runtime,
                    'url': url,
                    'date': gitdate
                })

    print("Generating '{}'".format(outpath))
    with open(outpath, 'w') as f:
        f.write("gData = [\n")
        for r in res:
            f.write("{:},\n".format(r))
        f.write("]")


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--history',
        type=str,
        default="../../sim/sim25_benchmark/history72/out_*_daint",
        help='Wildcard for out* files')
    parser.add_argument('--out',
                        type=str,
                        default="data.js",
                        help='Path to output')
    parser.add_argument(
        '--url',
        type=str,
        default="https://github.com/cselab/aphros/commits/",
        help='Prefix for url, to be appended by commit hash')
    args = parser.parse_args()
    gen(args.history, args.out, args.url)


if __name__ == "__main__":
    main()
