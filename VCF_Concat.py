#! usr/bin/env python3
# Author: Logan Wallace lwallac2@uoregon.edu, lwallac2@fredhutch.org
# Date: 12/21/2023

'''The purpose of this program is to concatenate variant call files from multiple individuals into a single bed file that can be loaded to a custom track on the ucsc genome browser.'''

# Module import
import re
import argparse
import os
import logging
import gzip
import pysam
import textwrap

# Set up some error logging features
logging.basicConfig(filename='VCF_Concat.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s', level=logging.DEBUG)

# Set some command line arguments
parser = argparse.ArgumentParser(
                    prog='VCF_Concat.py',
                    description='The purpose of this program is to concatenate variant call files from multiple individuals into a single bed file that can be loaded to a custom track on the ucsc genome browser.',
                    epilog='Good Luck!')

parser.add_argument("-d", "--directory", help = "The directory containing the VCF files which will be concatenated together.", type = str, default = "/fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation/")

parser.add_argument("-f", "--filenames", help = "A text file containing the filenames for the VCF files which will be concatenated together. One filename per line.", type = str, default = "/fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation/vcf_filenames.txt")

parser.add_argument("-m", "--file_mapping", help = "A text file containing a mapping of the sample names to our internal sample IDs (USIs). One sample name and USI per line, separated by a tab.", type = str, default = "/fh/fast/meshinchi_s/workingDir/scripts/lwallac2/Python/PacBio_LR_Initial_Investigation/file_mapping.txt")

parser.add_argument("-o", "--output", help = "The prefix for the output files. Will default to [variant_type].bed", type = str)

parser.add_argument("-c", "--collapse", help = "If this flag is set, the program will collapse the nearby variant calls of an identical SVTYPE into a single record in the output file. This is useful for visualizing structural variants on the UCSC genome browser. Defaults to True. To stop this collapsing, --collapse False", type = bool, default = True)

parser.add_argument("-w", "--window", help = "The window size for collapsing variant calls of the same SVTYPE. Defaults to 1000bp.", type = int, default = 10)

parser.add_argument("-r", "--records", help = "The number of records to keep in the LIFO list. Defaults to 5.", type = int, default = 5)

parser.add_argument("-g", "--gtf", help = "The name of the GTF file to use for the gene annotations. The file must be sorted, compressed and indexed using tabix. Defaults to gencode.v38.annotation.sorted.gtf.gz", type = str, default = "gencode.v38.annotation.sorted.gtf.gz")

args = parser.parse_args()

# Variable declaration
records = args.records
window = args.window
names_file = args.filenames 
directory = args.directory
SVTYPEs = ['BND', 'DEL', 'DUP', 'INS', 'INV']
BND_Interchrom_Color = ('55,126,184')
BND_Intrachrom_Color = ('228,26,28')
DUP_Color = ('255,127,0')
DEL_Color = ('152,78,163')
INS_Color = ('77,175,74')   
INV_Color = ('166,86,40')
HG_38_Chrom_List = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 
                    'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 
                    'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                    'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

# Class definitions
class LifoList:
    '''This class creates a list with a maximum size and a LIFO (last in, first out) structure that can be used to check records for collapsing. See argument - records'''
    def __init__(self, max_size=5):
        self.max_size = max_size
        self.list = []

    def add(self, item):
        # Add item at the beginning
        self.list.insert(0, item)
        # Ensure the list doesn't exceed the max size
        if len(self.list) > self.max_size:
            self.list.pop()

    def __str__(self):
        return str(self.list)


# Functions
def get_filenames(names_file):
    '''This function takes a text file containing the names of the VCF files to be concatenated together and returns a list of the filenames.'''
    with open(names_file, 'r') as f:
        filenames = f.readlines()
        filenames = [filename.strip() for filename in filenames]
    return filenames

def check_files(filenames):
    '''This function takes a list of filenames and checks to see if the files exist in the directory specified by the user. If the file does not exist, the program will exit with an error message.'''
    for filename in filenames:
        if os.path.isfile(filename):
            pass
        else:
            logging.error("File {} missing.".format(filename))
            exit()

def get_file_mapping(file_mapping):
    '''This function will take a filename as input and return a dictionary mapping the sample name to our internal sample ID (USI). This is an optional step that will be used if the user wants to convert the sample names to USIs.'''
    with open(file_mapping, 'r') as f:
        mapping = f.readlines()
        mapping = [line.strip().split('\t') for line in mapping]
        mapping = {line[0]:line[1] for line in mapping}
    return mapping

def get_record(line):
    # Split the line into a list, saving each field as needed 
    line = line.split('\t')
    CHROM = line[0]
    POS = line[1]
    ID = line[2]
    REF = line[3]
    ALT = line[4]
    QUAL = line[5]
    FILTER = line[6]
    INFO = line[7]
    FORMAT = line[8]
    # Pull the newline character off the end of SAMPLE field
    SAMPLE = line[9].strip('\n')
    # Get the variant type because this informs how we will parse the rest of the line
    variant_type = line[7].split(';')[0].replace('SVTYPE=', '')
    # Generate the output line for each type of variant
    if variant_type == 'BND':
        chr = CHROM
        chromStart = POS
        # Initialize the gene variables
        gene = None
        gene2 = None
        # chromEnd depends on if inter or intra chromosomal. 
        # Extract the chromosome IDs from the ID field by searching for the numbers that directly follow the 'chr' part of the ID string
        match = re.search(r'chr(.+):(\d+)-chr(.+):(\d+)', ID)
        if match:
            chr1, pos1, chr2, pos2 = match.groups()
            chr1 = 'chr' + chr1
            chr2 = 'chr' + chr2
            # Check to see if the chromosome IDs are in the list of chromosomes for hg38
            if chr1 not in HG_38_Chrom_List or chr2 not in HG_38_Chrom_List:
                missing = [chr for chr in [chr1, chr2] if chr not in HG_38_Chrom_List]
                logging.error("Chromosome ID not recognized: {}".format(missing))
                record = ""
                return record, variant_type
            # Check to see if the variant is inter or intra chromosomal
            # convert the positions to integers
            if chr1 == chr2:
                # We need to set a rule for intrachromosomal fusions that if the position 2 is less than position 1 we do not write out the record because it is effectively a duplicate of the record with the positions swapped. 
                if int(pos2) < int(pos1):
                    record = ""
                    return record, variant_type
                elif int(pos2) > int(pos1):
                    # If intra chromosomal, chromEnd is the second position
                    chromEnd = pos2
                    itemRgb = ('228,26,28')
            else:
                # If inter chromosomal, chromEnd is the first position plus 5 
                pos1 = int(pos1)
                chromEnd = pos1 + 5
                chromEnd = str(chromEnd)
                itemRgb = ('55,126,184')            
            # For name, we want to have the USI, the variant type, and the variant ID and allelic information. Later, we may want to incorporate more complex information into the field name, especially for BNDs. This might include the partner fusion gene, etc. 
            # Later is now
            # Find the partner fusion gene by querying the gtf file 
            # Find the gene for position 1 
            # Initialize the gene variable
            for rec in gtf.fetch(chr1, int(pos1), (int(pos1) + 1)):
                rec = rec.split('\t')
                # Find the gene_name from the 8th field 
                match = re.search(r'gene_name \"(.+?)\";', rec[8])
                if match:
                    gene = match.groups()[0]
            # Find the gene for position 2 
            for rec in gtf.fetch(chr2, int(pos2), (int(pos2) + 1)):
                rec = rec.split('\t')
                # Find the gene_name from the 8th field 
                match = re.search(r'gene_name \"(.+?)\";', rec[8])
                if match:
                    gene2 = match.groups()[0]
            # If the gtf.fetch found a gene name, add it to the ID field
            if gene:    
                ID = ID + '_' + gene
            if gene2:
                ID = ID + '_' + gene2
            # Separate the sample field to retrieve each component 
            sample = SAMPLE.split(':')
            GT = sample[0]
            AD = sample[1]
            REF = int(AD.split(',')[0]) 
            ALT = int(AD.split(',')[1])
            DP = sample[2]
            # Define a ratio of the REF and ALT alleles
            if int(DP) != 0:
                fusion_perc = int(AD.split(',')[1]) / int(DP)
            else:
                fusion_perc = 1
            # Limit the fusion percentage to 2 decimal places
            fusion_perc = str(round(fusion_perc, 2))
            name = f"{USI}_{variant_type}_{ID}_REF:{REF}_ALT:{ALT}_FusionPercent:{fusion_perc}_GT:{GT}"
            score = "0"
            strand = '.'
            thickStart = chromStart
            thickEnd = chromEnd          
            record = ('\t'.join([chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb]) + '\n')                                                  
        else:
            logging.error("Could not parse the ID field for the BND variant: {}".format(ID))
            record = ""
            return record, variant_type
    elif variant_type == 'DEL':
        chr = CHROM
        if chr not in HG_38_Chrom_List:
            logging.error("Chromosome ID not recognized: {}".format(chr))
            record = ""
            return record, variant_type
        chromStart = POS
        # Get the end position from the INFO field
        match = re.search(r'END=(\d+)', INFO)
        if match:
            end = match.groups()
            chromEnd = end[0]
            # Separate the sample field to retrieve each component 
            sample = SAMPLE.split(':')
            GT = sample[0]
            AD = sample[1]
            REF = int(AD.split(',')[0]) 
            ALT = int(AD.split(',')[1])
            DP = sample[2]
            # Define a ratio of the REF and ALT alleles
            if int(DP) != 0:
                fusion_perc = int(AD.split(',')[1]) / int(DP)
            else:
                fusion_perc = 1
            # Limit the fusion percentage to 2 decimal places
            fusion_perc = str(round(fusion_perc, 2))
            name = f"{USI}_{variant_type}_{ID}_REF:{REF}_ALT:{ALT}_FusionPercent:{fusion_perc}_GT:{GT}"
            score = "0"
            strand = '.'
            thickStart = chromStart
            thickEnd = chromEnd
            itemRgb = DEL_Color
            record = ('\t'.join([chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb]) + '\n')
        else:
            logging.error("Could not parse the INFO field for the DEL variant: {}".format(ID))
            record = ""
            return record, variant_type
    elif variant_type == 'DUP':
        chr = CHROM
        if chr not in HG_38_Chrom_List:
            logging.error("Chromosome ID not recognized: {}".format(chr))
            record = ""
            return record, variant_type        
        chromStart = POS
        # Get the end position from the INFO field
        match = re.search(r'END=(\d+);', INFO)
        if match:
            end = match.groups()
            chromEnd = end[0]
            # Separate the sample field to retrieve each component 
            sample = SAMPLE.split(':')
            GT = sample[0]
            AD = sample[1]
            REF = int(AD.split(',')[0]) 
            ALT = int(AD.split(',')[1])
            DP = sample[2]
            # Define a ratio of the REF and ALT alleles
            if int(DP) != 0:
                fusion_perc = int(AD.split(',')[1]) / int(DP)
            else:
                fusion_perc = 1
            # Limit the fusion percentage to 2 decimal places
            fusion_perc = str(round(fusion_perc, 2))
            name = f"{USI}_{variant_type}_{ID}_REF:{REF}_ALT:{ALT}_FusionPercent:{fusion_perc}_GT:{GT}"
            score = "0"
            strand = '.'
            thickStart = chromStart
            thickEnd = chromEnd
            itemRgb = DUP_Color
            record = ('\t'.join([chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb]) + '\n')
        else:
            logging.error("Could not parse the INFO field for the DUP variant: {}".format(ID))
            record = ""
            return record, variant_type
    elif variant_type == 'INS':
        chr = CHROM
        if chr not in HG_38_Chrom_List:
            logging.error("Chromosome ID not recognized: {}".format(chr))
            record = ""
            return record, variant_type
        chromStart = POS
        # Get the end position from the INFO field
        match = re.search(r'END=(\d+);', INFO)
        if match:
            end = match.groups()
            chromEnd = end[0]
            # Separate the sample field to retrieve each component 
            sample = SAMPLE.split(':')
            GT = sample[0]
            AD = sample[1]
            REF = int(AD.split(',')[0]) 
            ALT = int(AD.split(',')[1])
            DP = sample[2]
            # Define a ratio of the REF and ALT alleles
            if int(DP) != 0:
                fusion_perc = int(AD.split(',')[1]) / int(DP)
            else:
                fusion_perc = 1
            # Limit the fusion percentage to 2 decimal places
            fusion_perc = str(round(fusion_perc, 2))
            name = f"{USI}_{variant_type}_{ID}_REF:{REF}_ALT:{ALT}_FusionPercent:{fusion_perc}_GT:{GT}"
            score = "0"
            strand = '.'
            thickStart = chromStart
            thickEnd = chromEnd
            itemRgb = INS_Color
            record = ('\t'.join([chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb]) + '\n')
        else:
            logging.error("Could not parse the INFO field for the INS variant: {}".format(ID))
            exit()
    elif variant_type == 'INV':
        chr = CHROM
        chromStart = POS
        # Get the end position from the INFO field
        match = re.search(r'END=(\d+);', INFO)
        if match:
            end = match.groups()
            chromEnd = end[0]
            # Separate the sample field to retrieve each component 
            sample = SAMPLE.split(':')
            GT = sample[0]
            AD = sample[1]
            REF = int(AD.split(',')[0]) 
            ALT = int(AD.split(',')[1])
            DP = sample[2]
            # Define a ratio of the REF and ALT alleles
            if int(DP) != 0:
                fusion_perc = int(AD.split(',')[1]) / int(DP)
            else:
                fusion_perc = 1
            # Limit the fusion percentage to 2 decimal places
            fusion_perc = str(round(fusion_perc, 2))
            name = f"{USI}_{variant_type}_{ID}_REF:{REF}_ALT:{ALT}_FusionPercent:{fusion_perc}_GT:{GT}"
            score = "0"
            strand = '.'
            thickStart = chromStart
            thickEnd = chromEnd
            itemRgb = INV_Color
            record = ('\t'.join([chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb]) + '\n')
        else:
            logging.error("Could not parse the INFO field for the INV variant: {}".format(ID))
            record = ""
            return record, variant_type
    else:
        # If the variant type is not recognized, send a message to the logging file and keep going
        logging.error("Variant type not recognized: {}".format(variant_type))
        record = ""
        variant_type = ""
    return record, variant_type


# Retrieve the filenames from the names_file
filenames = get_filenames(names_file)

# Check to see if the files exist in the directory specified by the user
check_files(filenames)

# Get the file mapping if the user wants to convert the sample names to USIs
if args.file_mapping:
    mapping = get_file_mapping(args.file_mapping)

# Open a dictionary to store the files for each SVTYPE
files = {}

# Open the output files for writing, one for each SVTYPE
for var in SVTYPEs:
    # If the user specified an output prefix, the output file will be named output_prefix.SVTYPE.bed
    if args.output:
        output = args.output + '.' + var + '.bed'
    # Otherwise, the output file will be named SVTYPE.bed
    else:
        output = var + '.bed'
    # Open each output file for writing as the variable SVTYPE
    files[var] = open(output, 'w')
    # Write the header line to the output file
    files[var].write(f"track name={var}_VCF_Concat description=\"VCF Concatenation for the {var} variants\" visibility=2 itemRgb=\"On\"\n")

# Open the GTF file for annotation of the BND records using pysam
gtf = pysam.TabixFile(args.gtf)
# Initalize the LIFO list
recent_records = LifoList()
# Iterate through the VCFs
for file in filenames:  
    # Open the file for reading, using gzip.open if it's a .gz file
    open_func = gzip.open if file.endswith('.gz') else open
    with open_func(file, 'rt') as f:            
        # Iterate through the lines in the file
        for line in f:
            # Skip the header lines
            if line.startswith('##'):
                continue
            # Find the format line
            if line.startswith('#CHROM'):
                # Remove the '#' from the beginning of the line
                line = line.strip('#').rstrip()
                # Split the line into a list
                line = line.split('\t')
                # If the user wants to convert the sample names to USIs, do so
                if args.file_mapping:
                    line[9:] = [mapping[sample] for sample in line[9:]]
                # Save the format of the line as a list
                format = line
                # Remove the "['" and "']" from the beginning and end of the USI list
                USI = line[9]
            # Otherwise, the line is a variant record, parse the record
            else:
                # Get the records
                record, variant_type = get_record(line)
                # Set a boolean flag to check for duplicates inside our LIFO list loop
                duplicate = False
                # If the collapse flag is set, collapse the nearby variant calls of an identical SVTYPE into a single record in the output file.
                if args.collapse:
                    # If the record is empty, skip it
                    if record == "":
                        continue
                    # Otherwise, check the record against the other records in our list
                    else:
                        # If the list is empty, write the record to the output file
                        if len(recent_records.list) == 0:
                            # Add the record to the list of recent records
                            recent_records.add(record)
                        elif len(recent_records.list) > 0:
                            # Iterate through the records in the list
                            for recent in recent_records.list:
                                # If the record is the same as one of the recent records, skip it
                                if record == recent:
                                    continue
                                # Otherwise, check to see if the chromosome, USI and SVTYPE are the same and the distance start positions are within the window size
                                else:
                                    record_check = record
                                    # Split the record into a list
                                    record_check = record_check.split('\t')
                                    # Split the recent record into a list
                                    recent = recent.split('\t')
                                    # If the start chromosome is the same
                                    # If the start position is within the window size
                                    if record_check[0] == recent[0] and abs(int(record_check[1]) - int(recent[1])) < args.window:
                                        # Split the name field on the literal '.' character to check the USI and SVTYPE
                                        record_name = record_check[3].split('.')
                                        recent_name = recent[3].split('.')
                                        # If the USI and SVTYPE are the same
                                        if record_name[0] == recent_name[0]:
                                            # Set the duplicate flag to True
                                            duplicate = True
                                            # Break out of the loop
                                            break
                            # If the record is a duplicate, write it out to the logging file and not the output file
                            if duplicate == True:
                                logging.info("Duplicate record: {}".format(record))
                            # Otherwise, write the record to the output file
                            elif duplicate == False:
                                # Write the record to the output file for the appropriate variant type
                                files[variant_type].write(record)
                                # Add the record to the list of recent records
                                recent_records.add(record)
                # If the collapse flag is not set, write the record to the output file for the specific variant type
                else:
                    files[variant_type].write(record)

# Print a complete statement to the user
print("Concatenation completed. Output files written to the current directory.")

# Close the log file
logging.shutdown()