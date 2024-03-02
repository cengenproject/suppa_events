#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=merge
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --time=23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexis.weinreb@yale.edu


set -ue

echo "Starting $(date)"


tmp_dir="/vast/palmer/scratch/hammarlund/aw853/240229_salmon"

in_dir="/gpfs/gibbs/pi/hammarlund/CeNGEN/bulk/fastq/bulk"



mkdir -p $tmp_dir


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

# no duplicates
samples=($(printf "%s\n" "${samples[@]}" | sort -u))

samples=("${samples[@]:150:75}")

echo "Processing ${#samples[@]} samples"

echo " -------------------------------- "



echo
echo "##############################################"
echo "#########  Transfer and merge files  #########"
echo "##############################################"
echo
  

for samp in "${samples[@]}"
do

  if [ -e "$tmp_dir/${samp}_I1.fq.gz" ]
  then
  
    echo "Skipping $samp"
  
  else

    echo "-----  Sample $samp "
    
    cat $in_dir/*/Sample_$samp/${samp}_*_R1_00?.fastq.gz > $tmp_dir/${samp}_R1.fq.gz
    cat $in_dir/*/Sample_$samp/${samp}_*_R1_00?.fastq.gz > $tmp_dir/${samp}_R2.fq.gz
  
    zcat $in_dir/*/Sample_$samp/${samp}_*_R2_00?.fastq.gz \
      | awk 'NR%4 == 1 {print} NR%4 == 2 {print substr($0, 0, 8)} NR%4 == 3 {print} NR%4 == 0 {print substr($0,0,8)}' \
      | gzip -c \
      > $tmp_dir/${samp}_I1.fq.gz
  
  fi

done
  
echo
echo "##############################################"
echo "#########         UMI extract        #########"
echo "##############################################"
echo
  
for samp in "${samples[@]}"
do
  
  echo "-----  Sample $samp "
  
  umi_tools extract \\
    -I $tmp_dir/${samp}_I1.fq.gz \\
    -S $tmp_dir/$samp.out \\
    --log $tmp_dir/$samp.log \\
    --read2-in=$tmp_dir/${samp}_R1.fq.gz \\
    --read2-out=$tmp_dir/${samp}_umi_R1.fq.gz \\
    --extract-method=string \\
    --error-correct-cell \\
    --bc-pattern=NNNNNNNN \\
    >> $samp.umiR1.log
  
  
  umi_tools extract \\
    -I $tmp_dir/${samp}_I1.fq.gz \\
    --log $tmp_dir/$samp.log \\
    --read2-in=$tmp_dir/${samp}_R2.fq.gz \\
    --read2-out=$tmp_dir/${samp}_umi_R2.fq.gz \\
    --extract-method=string \\
    --error-correct-cell \\
    --bc-pattern=NNNNNNNN \\
    >> $samp.umiR2.log
  
done
  
  
  



echo "Finished $(date)"



