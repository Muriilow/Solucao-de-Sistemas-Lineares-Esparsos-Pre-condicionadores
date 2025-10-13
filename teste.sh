#!/bin/bash

PROG=cgSolver
CPU=3

DATA_DIR="Dados/"

mkdir -p ${DATA_DIR}

make purge
make

METRICA="FLOPS_DP L3CACHE ENERGY"
TEMPOS="${DATA_DIR}/Tempos.csv"
INPUTS=("SEM_1.in" "COM_1.in" "SEM_2.in" "COM_2.in" "SEM_3.in" "COM_3.in" "SEM_4.in" "COM_4.in" "SEM_5.in" "COM_5.in" "SEM_6.in" "COM_6.in")

for m in ${METRICA}
do
    
    LIKWID_CSV="${DATA_DIR}/${m}.csv"
        rm -f ${TEMPOS}

    for input_file in "${INPUTS[@]}"; do
        LIKWID_OUT="${DATA_DIR}/${m}_${input}.txt"
        echo -e "$input_file"

        output_file="${input_file%.*}.out"

        ./"$PROG" <  "$DATA_DIR""$input_file" > "$DATA_DIR""$output_file"
        #expected_file="${input_file%.*}"
        #diff "$DATA_DIR""$output_file" "$DATA_DIR""$expected_file"

        likwid-perfctr -O -C ${CPU} -g ${m} -o ${LIKWID_OUT} -m ./${PROG} ${n} >>${TEMPOS}
    done
done

