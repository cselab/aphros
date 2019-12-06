#!/usr/bin/env bash
cores_per_node=12

# Histogram metrics
SYNCHID=2
stats=(
    "hist_distrmesh.bin Run"
    "hist_distrmesh.bin RunKernels(inner)"
    "hist_distrmesh.bin RunKernels(halo)"
    "hist_synch${SYNCHID}.bin avail_halo"
    "hist_synch${SYNCHID}.bin avail_inner"
    "hist_synch${SYNCHID}.bin halo_waitall"
    "hist_synch${SYNCHID}.bin init_maps"
    "hist_synch${SYNCHID}.bin pack"
    "hist_synch${SYNCHID}.bin compress_face"
    "hist_synch${SYNCHID}.bin halo_decompress"
)

get_fname()
{
    fhist="$1"; shift
    stat="$1"; shift
    fname="${stat}_${fhist}"
    echo "${fname%.*}.dat"
}

fheaders()
{
    header="# kp\tnodes\tcores\tmean\tsdev\tsum"
    cwd="$1"
    for s in "${stats[@]}"; do
        IFS=' ' read -r -a file_stat <<< "${s}"
        fname=$(get_fname "${file_stat[@]}")
        echo -e ${header} > "${cwd}/$fname"
    done
}

stats()
{
    cwd="$1"
    root="$(pwd -P)"
    for s in "${stats[@]}"; do
        IFS=' ' read -r -a file_stat <<< "${s}"
        fhist="${file_stat[0]}"
        stat="${file_stat[1]}"
        fname="$(get_fname "${file_stat[@]}")"
        cd "${cwd}"
        for n in $(ls -d kp* | sort -t '=' -k2n); do
            kp=${n##*=}
            nodes=$((kp * kp * kp))
            cores=$((nodes * cores_per_node))
            if [[ ! -f "${n}/${fhist}" ]]; then
                continue
            fi
            ch.histbin -f "${n}/${fhist}" -rs -sn ${stat} > 000_tmp.dat
            mean="$(cat 000_tmp.dat | awk '/^\s*mean/{ print $3; exit }')"
            sdev="$(cat 000_tmp.dat | awk '/^\s*sdev/{ print $3; exit }')"
            sum="$(cat 000_tmp.dat | awk '/^\s*sum/{ print $3; exit }')"
            echo -e "${kp}\t${nodes}\t${cores}\t${mean}\t${sdev}\t${sum}" >> $fname
            rm 000_tmp.dat
        done
        cd "${root}"
    done
}

root="$(pwd -P)"
cwd='weak_scaling'
if [[ $# -eq 1 ]]; then
    cwd="${1}"
fi

for c in ${cwd}/omp*; do
    fheaders ${c}
    stats ${c}
done

# gnuplot scripts
cd ${cwd}
for s in "${stats[@]}"; do
    IFS=' ' read -r -a file_stat <<< "${s}"
    fname=$(get_fname "${file_stat[@]}")
    fgp="${fname%.*}.gp"
    fhist="${file_stat[0]}"
    stat="${file_stat[1]}"
    cat <<EOF > ${fgp}
#!/usr/bin/env gnuplot
set term x11 noenhanced
set key top right
set grid
set title '${fhist} -- ${stat}'
set xrange [0:11]
set xlabel 'kp'
set xtics 0, 1, 11
set ylabel 'seconds'
EOF
    chmod 755 "${fgp}"
    for c in omp*; do
        echo "replot '${c}/${fname}' using 1:6 w p ps 1 t '${c}'" >> "${fgp}"
    done
    sed -i 's/\(^replot.*$\)/\1, \\/' "${fgp}"
    sed -i '0,/replot/s/replot/plot/' "${fgp}"
    tac ${fgp}  | sed '0,/replot/s/, \\$//' > 000_tmp
    tac 000_tmp | sed 's/replot/ /' > "${fgp}"
    cat <<EOF >> "${fgp}"
pause(-1)
EOF
done
cd ${root}
