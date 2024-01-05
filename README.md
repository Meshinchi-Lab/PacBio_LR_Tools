# PacBio Long Read Sequencing Analysis Tools

Tools for the analysis of Pacific Biosciences long read sequencing data

### Overview 

Tools developed for the purpose of analyzing the long read WGS and RNAseq data resulting from the [Gabriella Miller Kids First X01 Project](https://commonfund.nih.gov/kidsfirst/2023X01projects#FY23_Meshinchi). Investigation of germline and somatic variants relevant to childhood AML, Down-Syndrome AML and TAM. 

Following sequencing, data was analyzed according to this [workflow](https://github.com/kids-first/kf-longreads-workflow/blob/main/docs/PACBIO_WORKFLOW_README.md) and is hosted on the Cavatica platform [here](https://cavatica.sbgenomics.com/u/kids-first-drc/sd-pet7q6f2/files/#q). 

## Tools

### Variant Call File Concatenation and BED File Generation - VariantVista

<p align="center">
  <img src="https://github.com/Meshinchi-Lab/PacBio_LR_Tools/blob/main/VariantVista.jpeg">
</p>

The application of various structural variant calling algorithms ([PBSV](https://github.com/PacificBiosciences/pbsv), [Sniffles](https://github.com/fritzsedlazeck/Sniffles), [longreadsv](https://support.sentieon.com/manual/)) to both the whole genome and RNA sequencing long read data resulted in a large amount of variant call files (VCFs). With a lack of visualization options, we decided to write a custom python script which will take a series of VCF files and generate BED files for each structural variant type that can be uploaded to a custom track in the UCSC Genome Browser and allow for the simultaneous visualization of all the structual variants at one time. 

VariantVista is intended to run on a collection of VCF files. 

#### Prerequisites

- Python `3.7.4`
- tabix `0.2.6`

#### Installation

To download this repository - 

```bash
git clone https://github.com/Meshinchi-Lab/PacBio_LR_Tools.git
```

#### VariantVista Usage 

Listing the command line arguments for VariantVista

```bash
python VariantVista.py --help
```

Running VariantVista.py on a collection of VCF files to produce a series of BED files able to be loaded to the UCSC Genome Browser as a custom track.
- -f vcf_filenames.txt in this repository for an example
- -m file_mapping.txt in this repository for an example

```bash
python VariantVista.py -d /directory/containing/VCF/files -f VCF_Filenames.txt -m file_mapping.txt -o example_run -c TRUE -w 10 -r 10 -g gencodev38.annotation.sorted.gtf.gz
```

## TODO

1. Add a flag to records already written out, that will notify the user that a near identical variant has been excluded from the BED file and the visualization.
2. Make the code more robust to different variant caller formats. 
3. Add argument to make the gene annotation with a GTF file optional. 



