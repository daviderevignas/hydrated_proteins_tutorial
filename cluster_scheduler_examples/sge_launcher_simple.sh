#!/bin/bash
for ((i=1; i<=60; i++))
do
   qsub -l mem_free=3G new_jobV.sh   
done
