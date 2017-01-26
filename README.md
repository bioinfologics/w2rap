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
* [QUAST] (http://quast.sourceforge.net/quast)

## w2rap steps using Saccharomyces cerevisiae dataset
### QC PE read files
1) Run FASTQC to check read qualities etc.

```
mkdir fastqc
fastqc -o fastqc scer_R1.fastq scer_R2.fastq
```
You can view the HTML reports generated in the fastqc directory.

2) Calculate read count and coverage;  
We have 3,648,316 PE reads (read length 100bp)  
So we have 3,648,316 * 100 * 2 = 729,663,200 bp coverage   
The [S. cerevisiae genome] (http://www.biology-pages.info/G/GenomeSizes.html) is ~12.5 Mb  
This equates to 729,663,200 / 12,495,682 = 58.4x genome coverage
 
3) Use KAT hist to generate a kmer histogram to estimate kmer coverage

```
kat hist -o scer_pe_hist -h 80 -t 8 -m 27 -H 100000000 scer_R?.fastq
```
![] (images/scer_pe_hist.png)

4) Use KAT comp to create a density plot comparing read 1 and read 2

```
kat comp -o scer_pe_R1vsR2 -n -t 8 -m 27 -H 100000000 -I 100000000 scer_R1.fastq scer_R2.fastq
```
![] (images/scer_pe_R1vsR2-main.mx.density.png)

5)  Download the S. cerevisiae [reference] (http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/), map reads and generate a SAM file. 

```
bwa index -p scer_ref -a bwtsw ref/S288C_reference_sequence_R64-2-1_20150113.fsa
bwa mem -SP -t 8 scer_ref scer_R?.fastq > pe2ref.sam
```

6) Generate an insert size histogram to check the insert size and shape of the distribution.

```
grep -v ‘@SQ' pe2ref.sam | grep -v '@PG' | awk -v binsize=20 '{if ($5>40) {if ($9>0) {print int($9/binsize)}else{print int($9/binsize*-1)}}}' | sort -n | uniq -c | awk -v binsize=20 '{print $2*binsize","$1}' > pe2ref.is
```
![](images/scer_pe_isize.png)
Insert size is 270bp

### Contigging
Use the w2rap-contigger to generate contigs from the PE reads.

```
mkdir contigs
w2rap-contigger/bin/w2rap-contigger -t 16 -m 200 -r scer_R1_san.fastq,scer_R2_san.fastq -o contigs -p scer_k200 -K 200
```
The contigs FASTA is generated in contigs/a.lines.fasta 

### Contig assessment
1) Check N50, total content etc.

```
abyss-fac contigs/a.lines.fasta
```
![](images/contigs_fac.png)

2) Use KAT comp to compare PE reads to contigs

```
kat comp -o scer_pe_v2_ctgs -t 8 -m 27 -H 100000000 -I 100000000 'scer_R?.fastq' contigs/a.lines.fasta
```
![](images/scer_pe_v2_ctgs-main.mx.spectra-cn.png)

3) Align genes, QUAST, BUSCO etc.

### LMP processing
Run FastQC to check read qualities etc.

Run Python script to remove Nextera adapters from LMP reads and any PE contamination.  

```  
lmp_processing <read_file_list> <ncpus>  
```

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
* Align genes, QUAST, BUSCO etc.

### Generate release
* Check for contamination
* Remove phiX and Illumina adapters
