# w2rap
WGS (Wheat) Robust Assembly Pipeline

## Essential tools
KAT  
BWA (or other short-read aligner)  
FLASH  
FASTX toolkit  
Nextclip  

## Other tools
FASTQC  
BUSCO

## w2rap steps
### QC PE read files
* Run FASTQC etc. to check read qualities
* Calculate read count and coverage
* Use KAT hist to generate a kmer histogram to estimate kmer coverage
* Use KAT comp to compare read 1 and read 2
* Map the reads to a reference and generate an insert size histogram to check the insert size and the tightness of the distribution

### Contigging
How to run using multiple libraries?

### Contig assessment
* Check N50, total content, gaps etc.
* Use KAT comp to compare PE reads to contigs
* Align genes, BUSCO etc.

### LMP processing
Python script to remove the Nextera adapter from LMP reads and remove PE contamination from LMP libraries.  
  
Run as: lmp_processing \<read\_file\_list\> \<ncpus\>  

read\_file\_list: a text file containing a list of LMP FASTQ files to process.  Files must be uncompressed and end in \_R1.fastq or \_R2.fastq.  
ncpus: the number of CPUs to use.

eg.  

```  
/path/to/LIB1_R1.fastq  
/path/to/LIB1_R2.fastq  
/path/to/LIB2_R1.fastq  
/path/to/LIB2_R2.fastq  
```

Processed LMP files will be written to the 'nextclip' directory.

### QC processed LMPs
* Calculate read count and coverage
* Use KAT comp to compare LMP reads to PE reads
* Map the reads to a reference and generate an insert size histogram to check the insert size and the tightness of the distribution

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
