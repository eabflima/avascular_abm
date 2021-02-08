#!/bin/bash
rm *# *~ merged*.pdf &> /dev/null
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Plotting details -----------------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
terminal="set terminal pdf font 'Times-New-Roman,18' lw 1"
ext=pdf
font=23
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
#--------------- Plotting every scenario ----------
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
for file in dead*.dat; do
    samples=$(echo ${file} | cut -d_ -f2 | cut -d. -f1)
    std=$(tail ${file} -n1 | awk '{printf $2}')
    name=evol_${samples}
    echo ${terminal} > figura.cmd
    echo "set output '${name}.${ext}'" >> figura.cmd
    echo "set border linewidth 3" >> figura.cmd
    echo "set ytics nomirror" >> figura.cmd
    echo "set key ins vert reverse top Left left font ',${font}'" >> figura.cmd
    echo "set ylabel 'Confluence (Dead)' font ',${font}'" >> figura.cmd
    echo "set y2label 'Confluence (Live)' font ',${font}'" >> figura.cmd
    echo "set xlabel 'Time (hours)' font ',${font}'" >> figura.cmd
    echo "set ytics font ',${font}'" >> figura.cmd
    echo "set xtics 12 font ',${font}'" >> figura.cmd
    echo "set y2tics font ',${font}'" >> figura.cmd
    let fcap=font-5
    echo "set title '${samples} samples' offset 0,-0.75 font ',${fcap}'" >> figura.cmd
    echo -n "plot [-2:98][] " >> figura.cmd
    if [ ${std} = 'nan' ]
    then
        echo -n "'${file}' u (\$0*3):1 axis x1y1 t'' with p pt 7 ps 0.3 lc 'red'," | tee -a figura.cmd &> /dev/null
        echo -n "'live_${samples}.dat' u (\$0*3):1 axis x1y2 t'' with p pt 7 ps 0.3 lc 'blue'," | tee -a figura.cmd &> /dev/null
    else
        echo -n "'${file}' u (\$0*3):1:2 axis x1y1 t'' with yerrorbars pt 7 ps 0.3 lc 'red'," | tee -a figura.cmd &> /dev/null
        echo -n "'live_${samples}.dat' u (\$0*3):1:2 axis x1y2 t'' with yerrorbars pt 7 ps 0.3 lc 'blue'," | tee -a figura.cmd &> /dev/null
    fi
    gnuplot "figura.cmd"
    pdfcrop ${name}.pdf ${name}t.pdf &> /dev/null
    mv ${name}t.pdf ${name}.pdf
done
pdftk evol_*.pdf cat output merged_samples.pdf
rm evol_*.pdf figura.cmd &> /dev/null
