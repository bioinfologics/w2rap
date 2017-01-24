# w2rap
WGS (Wheat) Robust Assembly Pipeline

## Software required
To run the pipeline you will need to install the following;  

* [KAT] (https://github.com/TGAC/KAT)  
* [BWA] (https://sourceforge.net/projects/bio-bwa/files/) (or other short-read aligner)  
* [FLASh] (https://ccb.jhu.edu/software/FLASH/)  
* [FASTX toolkit] (http://hannonlab.cshl.edu/fastx_toolkit/)  
* [Nextclip] (https://github.com/richardmleggett/nextclip/)  
* Something to calculate assembly stats (eg. [abyss-fac] (http://www.bcgsc.ca/platform/bioinfo/software/abyss))

Other tools are optional depending on how much QC and validation you want to perform on your reads and assembly.  We recommend;  

* [FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
* [BUSCO] (http://busco.ezlab.org/)

## w2rap steps
### QC PE read files
* Run FASTQC to check read qualities etc.
* Calculate read count and coverage
* Use KAT hist to generate a kmer histogram to estimate kmer coverage
* Use KAT comp to compare read 1 and read 2
* Map the reads to a reference and generate an insert size histogram to check the insert size and shape of the distribution

### Contigging

### Contig assessment
* Check N50, total content, gaps etc.
* Use KAT comp to compare PE reads to contigs
* Align genes, BUSCO etc.

### LMP processing
Python script to remove Nextera adapters from LMP reads and remove any PE contamination.  
  
Run as: lmp_processing \<read\_file\_list\> \<ncpus\>  

read\_file\_list: a text file containing a list of LMP FASTQ files to process.  Files must be uncompressed and end in \_R1.fastq or \_R2.fastq.  
eg.  

```  
/path/to/LIB1_R1.fastq  
/path/to/LIB1_R2.fastq  
/path/to/LIB2_R1.fastq  
/path/to/LIB2_R2.fastq  
```

ncpus: the number of CPUs to use.

Processed LMP files will be written to the 'nextclip' directory. Read counts before and after trimming are written to the log file.

### QC processed LMPs
* Calculate fragment coverage from trimmed read count
* Use KAT comp to check for LMP representation issues by comparing LMP reads to PE reads
* Map the reads to a reference and generate an insert size histogram to check the insert size and the shape of the distribution

### Scaffolding
* s_prepare
* s_map (use PE and LMP reads)
* s_scaffold
* N-remapping script

### Scaffold validation
* Check N50, total content, gaps etc.
* Use KAT comp to compare PE reads to scaffolds
* Align genes, BUSCO etc.

### Generate release
* Check for contamination
* Remove phiX and Illumina adapters
