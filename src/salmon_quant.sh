#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=salm_quant
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


set -ue

echo "Starting $(date)"


tmp_dir="/vast/palmer/scratch/hammarlund/aw853/240229_salmon"

in_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/fastq/bulk"


sample_paths=($(echo $in_dir/*/Sample_*))
samples=()

for samp_path in "${sample_paths[@]}"
do
  if [[ "$samp_path" =~ /Sample_([A-Z0-9ef]{1,4}r[0-9]{1,4})$ ]];  then
    samples+=(${BASH_REMATCH[1]})
  else
    echo "Did not match: $samp_path "
  fi
done


echo "Processing ${#samples[@]} samples"

echo " -------------------------------- "



  
  
  

echo
echo "##############################################"
echo "#########       Salmon mapping       #########"
echo "##############################################"
echo
  

for samp in "${samples[@]}"
do
  
  
  
done


echo
echo "##############################################"
echo "#########       Salmon mapping       #########"
echo "##############################################"
echo
  

for samp in "${samples[@]}"
do

umi_tools dedup \\
                  -I aligned.bam \\
                  --paired \\
                  --output-stats=deduplicated \\
                  -S deduplicated.bam 




done



echo "Finished $(date)"



