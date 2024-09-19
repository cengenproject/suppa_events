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


events_file="data/events/${WS}_all_strict.ioe"
tx_tpm="data/240913_str_q_tx_TPM_filt.tsv"

out_dir="data/240918_psiPerEvent"

mkdir -p $out_dir


# Note we gathered all events in single file with:
# head -1 WS289_A3_strict.ioe > ../WS289_all_strict.ioe
# cat *.ioe | sed '/^seqname/d' >> ../WS289_all_strict.ioe


suppa.py psiPerEvent \
    --mode INFO \
    --ioe-file $events_file \
    --expression-file $tx_tpm \
    --total-filter 5 \
    --output-file $out_dir


echo "End $(date)"




