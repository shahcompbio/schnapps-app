#!/bin/bash
#BSUB -J schnapps
#BSUB -n 5
#BSUB -R rusage[mem=20]
#BSUB -W 4:00
#BSUB -eo %J.stderr

module load singularity/3.6.2

singularity exec --bind=/juno \
    /juno/work/shah/users/william1/projects/schnapps-app/schnapps_latest.sif \
    Rscript /juno/work/shah/users/william1/projects/schnapps-app/run-schnapps.R \
    --hmmcopyqc /juno/work/shah/isabl_data_lake/analyses/88/13/8813/results/A90554A_metrics.csv.gz /juno/work/shah/isabl_data_lake/analyses/88/14/8814/results/A90554B_metrics.csv.gz /juno/work/shah/isabl_data_lake/analyses/88/15/8815/results/A96213A_metrics.csv.gz \
    --hmmcopyreads /juno/work/shah/isabl_data_lake/analyses/86/78/8678/results/A96213A_reads.csv.gz /juno/work/shah/isabl_data_lake/analyses/86/79/8679/results/A90554B_reads.csv.gz /juno/work/shah/isabl_data_lake/analyses/86/80/8680/results/A90554A_reads.csv.gz \
    --allelecounts /juno/work/shah/isabl_data_lake/analyses/23/89/12389/results/allele_counts.csv.gz /juno/work/shah/isabl_data_lake/analyses/23/90/12390/results/allele_counts.csv.gz /juno/work/shah/isabl_data_lake/analyses/23/91/12391/results/allele_counts.csv.gz \
    --ncores 5 \
    --qcplot /juno/work/shah/users/william1/projects/schnapps-app/results/qcplot.png \
    --heatmap /juno/work/shah/users/william1/projects/schnapps-app/results/heatmap.png \
    --csvfile /juno/work/shah/users/william1/projects/schnapps-app/results/hscn.csv.gz \
    --Rdatafile /juno/work/shah/users/william1/projects/schnapps-app/results/schnapps.Rdata
