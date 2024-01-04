for file in VEP_Output/*.vcf.gz; do
  bs_name=$(zcat "$file" | grep '^#CHROM' | awk '{print $NF}')
  echo "$(basename $file): $bs_name" >> names.txt
done