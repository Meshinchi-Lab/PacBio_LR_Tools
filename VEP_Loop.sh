# Iterate over the array of .vcf.gz files
for vcf_file in VCF_Concat/*.vcf.gz; 
do
    sbatch VEP_Command.sh "$vcf_file"
done