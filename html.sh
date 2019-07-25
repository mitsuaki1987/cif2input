#!/bin/bash

for dir in *
do
    if [ -d ${dir} ]; then
        name=`echo ${dir} | sed -e 's/2s/2\//g' -e 's/3s/3\//g' -e 's/4s/4\//g' -e 's/6s/6\//g'`
        dos=`awk 'NR==2{printf "%f", $2}' ${dir}/${dir}.pdos_tot`
        printf "%s (%.3f [states/eV/u.c./spin]) </br>\n" ${name} ${dos}
        echo '<table border="0"><tr align="center">'
        for p in `find ${dir} -name "*.frmsf" |awk -F/ '{print $NF}'`
        do
            #printf '<td><a href="./%s" download="%s"><img src="./%s" height="200"></a></td>\n' \
            printf '<td><a href="./%s"><img src="./%s" height="200"></a></td>\n' \
                   ${dir}${p} ${dir}${p%frmsf}png
        done
        echo '</tr><tr align="center">'
        for p in ${dir}/${dir}.pdos_*[1-9]?
        do
            proj=${p#${dir}/${dir}.pdos_}
            pdos=`awk 'NR==2{printf "%f", $2}' ${p}`
            pdos=`echo ${pdos}/${dos}*100|bc -l`
            printf "<td>%s (%.3f %%)</td>\n" ${proj} ${pdos}
        done
        echo "</tr></table><hr>"
    fi
done
