# Logan Wallace
# 12.13.2023
# This loops through all of the VEP output files and compresses them with bgzip and indexes them with tabix for input to bcftools merge

for file in VEP_Output/*.vcf; do
  bgzip "$file"
  tabix -p vcf "$file.gz"
done