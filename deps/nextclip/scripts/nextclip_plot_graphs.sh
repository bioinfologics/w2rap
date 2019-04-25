#!/bin/bash

# Script:  nextclip_plot_graphs.sh
# Purpose: Plot insert length and read length graphs for NextClip
# Author:  Richard Leggett
# Contact: richard.leggett@tgac.ac.uk
#          http://www.tgac.ac.uk/richard-leggett/

libdir=$1
lib=$2
log=$3
scriptdir=$4

readsdir=${libdir}/reads
graphsdir=${libdir}/graphs
analysisdir=${libdir}/analysis

inputfile=${readsdir}/${lib}_duplicates.txt
outputfile=${graphsdir}/${lib}_duplicates.pdf
echo "" | tee -a ${log}
echo "Plotting duplication rate graph" | tee -a ${log}
echo " Input: ${inputfile}" | tee -a ${log}
echo "Output: ${outputfile}" | tee -a ${log}

if [ -f ${inputfile} ] ; then
    Rscript ${scriptdir}/nextclip_plot_duplication.R ${inputfile} ${outputfile}
else
echo "ERROR: nextclip_plot_graphs.sh: Can't find ${inputfile}" | tee -a ${log}
fi

for read in 1 2
do
    inputfile=${readsdir}/${lib}_R${read}_gc.txt
    outputfile=${graphsdir}/${lib}_R${read}_gc.pdf
    echo "Plotting GC content for read ${read}" | tee -a ${log}
    echo " Input: ${inputfile}" | tee -a ${log}
    echo "Output: ${outputfile}" | tee -a ${log}

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_gc.R ${inputfile} ${outputfile}
    else
        echo "ERROR: Can't find ${inputfile}" | tee -a ${log}
    fi
done

for type in A B C D
do
    inputfile=${readsdir}/${lib}_${type}_pair_hist.txt
    outputfile=${graphsdir}/${lib}_cumulative_pairs_${type}.pdf
    echo "Plotting cumulative pair lengths for category ${type}" | tee -a ${log}
    echo " Input: ${inputfile}" | tee -a ${log}
    echo "Output: ${outputfile}" | tee -a ${log}

    if [ -f ${inputfile} ] ; then
        Rscript ${scriptdir}/nextclip_plot_pair_lengths.R ${inputfile} ${outputfile}
    else
        echo "ERROR: nextclip_plot_graphs.sh: Can't find ${inputfile}" | tee -a ${log}
    fi

    for read in 1 2
    do
        inputfile=${readsdir}/${lib}_${type}_R${read}_hist.txt
        outputfile=${graphsdir}/${lib}_lengths_${type}_R${read}.pdf
        echo "Plotting lengths for category ${type} read ${read}" | tee -a ${log}
        echo " Input: ${inputfile}" | tee -a ${log}
        echo "Output: ${outputfile}" | tee -a ${log}

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_lengths.R ${inputfile} ${outputfile}
        else
            echo "ERROR: nextclip_plot_graphs.sh: Can't find ${inputfile}" | tee -a ${log}
        fi
    done

    for kind in mp pe tandem
    do
        inputfile=${analysisdir}/${lib}_${type}_${kind}.txt
        outputfile=${graphsdir}/${lib}_${kind}_${type}.pdf
        echo "Plotting insert size lengths for category ${type} kind ${kind}" | tee -a ${log}
        echo " Input: ${inputfile}" | tee -a ${log}
        echo "Output: ${outputfile}" | tee -a ${log}

        if [ -f ${inputfile} ] ; then
            Rscript ${scriptdir}/nextclip_plot_inserts.R ${inputfile} ${outputfile}
        else
            echo "ERROR: nextclip_plot_graphs.sh: Can't find ${inputfile}" | tee -a ${log}
        fi
    done
done
