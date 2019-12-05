#!/usr/bin/env bash
BS=16
LAYERS=4
NTHREADS=12

# if not using numa nodes on Euler but whole node
# myhost=$(hostname)
# if [[ ${{myhost:0:3}} == "eu-" ]]; then
#     # euler
#     NTHREADS=24
# fi

# full range (possibly too many)
kp=(1 2 3 4 5 6 7 8 9 10)

cases=(
    "omp=1_cubismnc=0_compress=0_layers=${LAYERS}_kb=1"           # pre-pasc

    "omp=1_cubismnc=1_compress=0_layers=${LAYERS}_kb=1"           # without field copying
    "omp=1_cubismnc=1_compress=0_layers=${LAYERS}_kb=2"
    "omp=1_cubismnc=1_compress=0_layers=${LAYERS}_kb=3"

    "omp=${NTHREADS}_cubismnc=1_compress=0_layers=${LAYERS}_kb=1" # without field copying + OpenMP
    "omp=${NTHREADS}_cubismnc=1_compress=0_layers=${LAYERS}_kb=2"
    "omp=${NTHREADS}_cubismnc=1_compress=0_layers=${LAYERS}_kb=3"

    "omp=1_cubismnc=1_compress=1_layers=${LAYERS}_kb=1"           # without field copying + compression
    "omp=1_cubismnc=1_compress=1_layers=${LAYERS}_kb=2"
    "omp=1_cubismnc=1_compress=1_layers=${LAYERS}_kb=3"

    "omp=${NTHREADS}_cubismnc=1_compress=1_layers=${LAYERS}_kb=1" # without field copying + OpenMP + compression
    "omp=${NTHREADS}_cubismnc=1_compress=1_layers=${LAYERS}_kb=2"
    "omp=${NTHREADS}_cubismnc=1_compress=1_layers=${LAYERS}_kb=3"
)

root="$(pwd -P)/weak_scaling"
base='base'
mkdir -p ${root}/${base}
cp -t ${root}/${base} add.conf Makefile pargen sim.makefile std.conf
cd ${root}

# enable histograms for scaling analysis
echo "set int histogram 1" >> ${base}/add.conf

for c in "${cases[@]}"; do
    cdir=${root}/${c}
    mkdir -p ${cdir}
    cd ${cdir}
    for k in "${kp[@]}"; do
        nodes="kp=${k}"
        if [[ -d ${nodes} ]]; then
            echo "Skip existing: ${c}/${nodes}"
            continue
        fi
        cp -r ${root}/${base} ${nodes}
        cd ${nodes}
        IFS='_' read -r -a opt <<< "${c}_${nodes}"
        make "${opt[@]}" "${@}"
        cd ${cdir}
    done
    cd ${root}
done
