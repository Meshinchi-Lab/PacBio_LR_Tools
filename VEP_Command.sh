#!/bin/bash
#SBATCH --job-name=VEP_WALLACE
#SBATCH --output=./slurmout/result-%j.out
#SBATCH --error=./slurmout/result-%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00

# Load modules
ml VEP/103.1-GCC-10.2.0

# Relevant directory
cd /fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation

#vcf_file=$(basename "$1" .vcf.gz)

# Run the VEP command for each of the input VEP files
#vep -i "VCF_Concat/${vcf_file}.vcf.gz" -o "/fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation/VEP_Output/${vcf_file}.pbsv.vep.annotation.vcf" --vcf --cache --gene_phenotype --symbol --dir_cache ensembl-vep --cache_version 110

vcf_file="882c63fc-a7be-440c-bacd-ab4710cfbfa6.pbsv.vcf.gz"

vep -i VCF_Concat/${vcf_file} -o /fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation/VEP_Output/${vcf_file}.pbsv.vep.annotation.vcf --vcf --cache --gene_phenotype --symbol --dir_cache ensembl-vep --cache_version 110