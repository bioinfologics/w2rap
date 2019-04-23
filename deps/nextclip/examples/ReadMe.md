#Examples#

This directory contains some example data for testing NextClip. This consists of:
1. Inside the reads directory: 125,000 pairs of reads from a (not especially high quality) Streptomyces coelicolor Nextera LMP run.
2. Inside the reference directory: the S. coelicolor reference, from https://www.ebi.ac.uk/ena/data/view/AL645882.
3. Inside the configure directory, LIB468.txt is a parameters file for passing to the NextClip pipeline.

Please be aware that the data is a much reduced subset of a full library, so results may have a few quirks (for example, only one PE alignment in range for category A resulting in an unusual graph). This is just intended as an illustration of how to use the tool.

#Running NextClip standalone#

First compile NextClip as detailed in the manual. To run the tool, change into the examples directory:

cd examples

Then type:

../bin/nextclip -i reads/LIB1468part_ATCACG_L001_R1_001.fastq -j reads/LIB1468part_ATCACG_L001_R2_001.fastq -n 125000 -o output

#Running the NextClip analysis pipeline#

To run the NextClip pipeline, ensure you have BWA, R and LaTeX installed and available and that you have carried out all the other installation steps detailed in the manual.

First, index the reference:

cd references

bwa index -a bwtsw AL645882.fasta

../../scripts/nextclip_index_reference.pl AL645882.fasta

Then to run the analysis pipeline, type:

cd ..

../scripts/nextclip_lmp_analysis.pl -config configure/LIB1468.txt -scheduler none -bwathreads 1

If all works successfully, a new directory called LIB1468 should appear and inside that will be a latex directory containing a PDF file.

Running the pipeline on a 2013-spec MacBook Pro took approximately 85 seconds.
