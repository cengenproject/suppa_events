#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=salmon_ind
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu



set -ue
echo "Start $(date)"



WS="289"


refpath="/gpfs/ycga/project/ysm/hammarlund/aw853/references/WS"$WS
txome_tmp_path="data/intermediates/2024-02-29_txome_with_decoys.fa.gz"



# Get references
module load R


rcmd_genome='cat(wbData::wb_get_genome_path('$WS', Sys.getenv("refpath")))'
rcmd_transcriptome='cat(wbData::wb_get_transcriptome_path('$WS', Sys.getenv("refpath")))'

genome_path=$( env refpath=$refpath Rscript -e "$rcmd_genome" )
transcriptome_path=$( env refpath=$refpath Rscript -e "$rcmd_transcriptome" )



# Prepare decoys

grep "^>" \
  <( gunzip -c $genome_path ) \
  | cut -d " " -f 1 \
  > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

cat $transcriptome_path $genome_path > $txome_tmp_path



echo "Indexing..."


# Salmon index

module switch R Salmon

salmon index \
  --transcripts $txome_tmp_path \
  --index 240229_index \
  --decoys decoys.txt \
  --threads $SLURM_CPUS_PER_TASK




echo
echo "Done $(date)"


