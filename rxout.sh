#!/bin/bash

nat=$(grep "number of atoms/cell" rx.out|awk '{print $5}' | tail -n 1)

sed -n -e '1,/CELL_PARAMETERS/p' "${1}" > temp
# shellcheck disable=SC2129
grep -A 3 CELL_PARAMETERS rx.out | tail -n 3 >> temp
awk '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS/' "${1}" >> temp
grep -A "${nat}" ATOMIC_POSITIONS rx.out |tail -n "${nat}" >> temp
sed -n -e '/K_POINTS/,$p' "${1}" >> temp

mv temp "${1}"
