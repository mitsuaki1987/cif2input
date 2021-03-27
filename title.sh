#!/bin/bash

authtoken=

for colcode in `cat code.dat`
do
    idnum=`curl -X GET "https://icsd.fiz-karlsruhe.de/ws/search/expert?query=COLLECTIONCODE%3A${colcode}" \
                 -H "accept: application/json" \
                      -H "ICSD-Auth-Token: ${authtoken}" | \
                          jq  '.idnums'[0] | sed -e "s/\"//g"`
    echo ${idnum}

    curl -X GET "https://icsd.fiz-karlsruhe.de/ws/csv?idnum=${idnum}&windowsclient=false&listSelection=Title" \
         -H "accept: application/csv" \
         -H "ICSD-Auth-Token: ${authtoken}" | awk 'NR==2{print '${colcode}', $0}' >> title.txt
done
