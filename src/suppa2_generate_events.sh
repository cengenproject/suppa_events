#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=suppa_genEvents
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu

set -ue

echo "Start $(date)"


module load miniconda
conda activate SUPPA2


WS="WS289"
ref_dir="/gpfs/ycga/project/ysm/hammarlund/aw853/references/$WS"
ref_gtf=$ref_dir"/c_elegans.PRJNA13758.${WS}.canonical_geneset.gtf"

out_dir="data/events/240229_events"

mkdir -p $out_dir


# SUPPA2 generate all local (ioe) events

suppa.py generateEvents \
    --mode INFO \
    --event-type SE SS MX RI FL \
    --format ioe \
    --boundary S \
    --input-file $ref_gtf \
    --output-file $out_dir/${WS}


echo "End $(date)"

