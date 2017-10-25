#!/bin/sh
#$ -cwd
#$ -l mf=30G,h_vmem=40G,h_fsize=2G
## -m e -M dengdetian0603@gmail.com

## -t 1-5

Rscript ./TestRun.R 3 200 1000 50
