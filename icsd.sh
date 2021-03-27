#!/bin/bash -x

# $1 = Login ID
# $2 = Password

for i in x*
do
    token=`curl -i -X POST "https://icsd.fiz-karlsruhe.de/ws/auth/login" -H "accept: text/plain" -H "Content-Type: application/x-www-form-urlencoded" -d "loginid=$1&password=$2" | awk '$1=="ICSD-Auth-Token:"{print $2}'`
    echo $i $token
    curl -i -X GET "https://icsd.fiz-karlsruhe.de/ws/cif/multiple?`awk 'NR==1{printf "idnum=%d", $1} NR>1{printf "&idnum=%d", $1}' ${i}`&windowsclient=false&filetype=zip" -H "accept: application/cif" -H  "ICSD-Auth-Token: ${token}" --output ${i}.zip
    echo ""
    curl -i -X GET "https://icsd.fiz-karlsruhe.de/ws/auth/logout" -H  "accept: text/plain" -H "ICSD-Auth-Token: ${token}"
    echo ""
done
