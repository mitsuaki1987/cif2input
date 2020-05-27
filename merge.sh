#!/bin/bash

paste ${1} ${2} | awk '{printf "%f ", $2;for(i=4;i<=NF;i++){printf "%f ", $(i)};print ""}'

