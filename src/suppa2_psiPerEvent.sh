#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=suppa_psiEvent
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


events_dir="data/events/240229_events"
tx_tpm="data/231208_str_q_tx_TPM.tsv"

out_dir="data/240301_psiPerEvent"

mkdir -p $out_dir


# Note, we used the command from the wiki to gather all events in single file:
#awk '
#    FNR==1 && NR!=1 { while (/^<header>/) getline; }
#    1 {print}
#' *.ioe > WS289_all_strict.ioe
# see https://github.com/comprna/SUPPA/wiki/SUPPA2-tutorial



suppa.py psiPerEvent \
    --mode INFO \
    --ioe-file $events_dir/${WS}_all_strict.ioe \
    --expression-file $tx_tpm \
    --output-file $out_dir


echo "End $(date)"




