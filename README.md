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
* [bioawk](https://github.com/lh3/bioawk)
* Python with Biopython and Matplotlib installed [python2] (https://www.python.org/downloads/release/python-2711/)

Other tools are optional depending on how much QC and validation you want to perform on your reads and assembly.  We recommend;  

* [FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
* [BUSCO] (http://busco.ezlab.org/)
* [QUAST] (http://quast.sourceforge.net/quast)

This tutorial assumes that you are using a Linux machine. If you do not have access to a Linux machine, you will need to find equivalent tools which run on your operating system to complete some of the steps.

## w2rap steps using Saccharomyces cerevisiae dataset
### 1) QC PE read files
a) Run FASTQC to check read metrics.

```
mkdir fastqc
fastqc -o fastqc scer_pe_R1.fastq scer_pe_R2.fastq
```
 [FastQC] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  generates an HTML report in the fastqc directory. The report will give some indications about the quality of your reads, and the success of your sequencing run. Consult the documentation on the website, and the example reports from good and bad Illumina runs for further information. You should calculate the read coverage using the read count. 

![] (images/fastqc.png)

FastQC shows we have 2,519,142 PE reads of length 300bp providing 2,519,142 * 300 * 2 = 1,511,485,200 bp coverage   
The [S. cerevisiae genome] (http://www.biology-pages.info/G/GenomeSizes.html) is ~12.5 Mb which means we have 1,511,485,200 / 12,495,682 = 121x genome coverage
 
b) Use KAT hist to generate a kmer histogram to estimate kmer coverage. The histogram shows us how often kmers appear in reads 1 and 2, 

```
kat hist -o scer_pe_hist -h 80 -t 8 -m 27 -H 100000000 scer_pe_R?.fastq
```
<img src="images/scer_pe_hist.png"  width="500" height="400">

We can see that as the frequency approaches zero, the number of distinct kmers increases significantly, these kmers are from erroneous reads. There is a relatively symmetrical distribution, centered at approximately 80, with a reasonably small variance. Hence, the estimated kmer coverage is equal to 80.  

c)  To enable a more detailed assessment of the quality of the reads, download the S. cerevisiae S288C [reference] (http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/), map the reads to it, and generate a SAM file. 

```
bwa index -p scer_ref -a bwtsw ref/S288C_reference_sequence_R64-2-1_20150113.fsa
bwa mem -SP -t 8 scer_ref scer_pe_R?.fastq > pe2ref.sam
```

By checking that a reasonable percentage of your reads map to the reference, you can be conifdent that your sequencing data is of a reasonable quality. You can do this by entering the following command:

```
samtools view -F 4  pe2ref.sam | cut -f 1 | sort | uniq | wc -l
```

d) From the SAM file generated above, we can obtain the raw data needed to draw an insert size histogram: 

```
grep -v ‘@SQ' pe2ref.sam | grep -v '@PG' | awk -v binsize=20 '{if ($5>40) {if ($9>0) {print int($9/binsize)}else{print int($9/binsize*-1)}}}' | sort -n | uniq -c | awk -v binsize=20 '{print $2*binsize","$1}' > pe2ref.is
```

Then, simply use your favourite plotting tool to check the insert size and the shape of the distribution.

<img src="images/yeast_pe_insert_size.png" width="550" height="400">

We can see that the insert sizes are roughly symmetrically distributed around 500. The distribution is quite wide, so a lot of pairs will have an insert size which varies quite far from the average, but we should be able to obtain a reasonable assembly from these reads.


### 2) Contigging

Use the w2rap-contigger to generate contigs from the PE reads. The current version of the w2rap contigger runs in 8 steps: 


Step # | Description | Outputs
:---|---|---
1 | Read loading | binary-formatted reads
2 | Kmer counting | raw kmer data
3 | Build small k (k=60) graph from reads | small k graph, read paths
4 | Build large K graph from small k graph and reads | large K graph, read paths
5 | Clean large K graph | large K cleaned graph, read paths
6 | Local assemblies on the large K graph "gaps" | large K completed graph, read paths
7 | Graph simplification and PathFinder | large K simplified graph, read paths, raw/contig-lines GFA and fasta
8 | PE-scale scaffolding across gaps in the large K graph | large K simplified graph with jumps, read paths, raw/lines GFA and fasta

By default the contigger will run each of these steps in order, not dumping unnecessary intermediate files. Each step can be run individually, by specifying the `--from_step ` and `--to_step`. If you specify the `--to_step`, the contigger will automatically dump the output files from the specified step. To be able to run from any intermediate step, the preceeding steps must have been run with the `--dump_all` flag set. 

You need to create a new directory for the intermediate and output files. To run from start to finish with default assembly parameters, run: 

```
mkdir contigs
w2rap-contigger/bin/w2rap-contigger -t 16 -m 200 -r scer_pe_R1.fastq,scer_pe_R2.fastq -o contigs -p scer_k200 
```
The contigs FASTA is generated in contigs/a.lines.fasta 


The number of times a kmer must appear in the reads to be included in the small k graph can be controlled with the `--min_freq` parameter:

```
w2rap-contigger/bin/w2rap-contigger -t 16 -m 200 -r scer_pe_R1.fastq,scer_pe_R2.fastq -o contigs -p scer_k200 --min_freq 20
```

Ideally, `--min_freq` should be selected to remove most error kmers, and retain most kmers which are genuinely present in the genome of interest. This value can be determined with the help of the kmer histogram from step b) of the QC. 

In the above examples we use the default kmer length of 200 but you may want to generate assemblies using different kmer lengths and assess each one. We can vary the value of k used to build the large k graph with the `-K` option, like so:

```
w2rap-contigger/bin/w2rap-contigger -t 16 -m 200 -r scer_pe_R1.fastq,scer_pe_R2.fastq -o contigs -p scer_k200 -K 220 --from_step 3
```

More detail about these options, and descriptions of the other options, can be found in the full w2rap paper, or by running:

```
w2rap-contigger/bin/w2rap-contigger --help
``` 

### 3) Contig assessment
a) Check contiguity stats.

```
abyss-fac contigs/a.lines.fasta
```
![](images/contigs_fac.png)

We are assembling 11.78e6 bps of our 12Mb genome in contigs longer than 500bp.  The contig-N50 is 173 Kb. The expected number of contigs and N50 will vary significantly between genomes, in particular more complex and repetitive genomes may be more fragmented and hence have a lower N50. 

b) Use KAT comp to generate a spectra-cn to compare PE reads to contigs

```
kat comp -o scer_pe_v2_ctgs -t 8 -m 27 -H 100000000 -I 100000000 'scer_pe_R?.fastq' contigs/a.lines.fasta
```

<img src="images/pe_vs_contigs_k27-main.mx.spectra-cn.png"  width="500" height="400">

This spectra shows we are assembling almost all the content from the reads correctly with no evidence of missassembly.  There is some evidence of reads from the error distribution appearing in the assembly (see the [KAT documentation](https://kat.readthedocs.io/en/latest/) for more details on how to interpret KAT plots).

c) Assess assembly accuracy and completeness using QUAST and aligning BUSCO genes.

```
mkdir quast
python /path/to/quast.py -o ./quast -R ref/S288C_reference_sequence_R64-2-1_20150113.fsa -t 8 -f ../tutorial_runthrough/contigs/a.lines.fasta
```
When a reference is provided, QUAST generates a report containing useful statistics including an estimation of missassemblies:


Genome statistics	 | a.lines
-------------------- |---------------
Genome fraction (%)			  |	91.919
Duplication ratio			  |	1.033
Largest alignment			  |	47316
Total aligned length		  |	11535940
NGA50							  |	9669
LGA50							  |	381
Misassemblies					  |
misassemblies					  |28
Misassembled contigs length  |	826088
Mismatches					  |
mismatches per 100 kbp		  |3.19
indels per 100 kbp			  |1.56
N's per 100 kbp				  |	55.41
Statistics without reference |	
contigs						  | 	1723
Largest contig				  |	76827
Total length					  |	11732078
Total length (>= 1000 bp)	  | 11480625
Total length (>= 10000 bp)	  | 5930056
Total length (>= 50000 bp)	  | 282022


Run BUSCO like so:

```
python /path/to/busco2/BUSCO.py -o busco_pe -in contigs/a.line.fasta -l ~/busco_data/eukaryota -m genome -f
```

The proportion of BUSCOs present is assumed to be similar to the proportion of all genes present, so the summary table enables us to estimate how well the assembly captures the genetic content of the genome:

	Count		|       Type    
------------ | -----------------------------------
        407  |   Complete BUSCOs
        379  |   Complete and single-copy BUSCOs
        28  |   Complete and duplicated BUSCOs
        13   |   Fragmented BUSCOs
        9    |   Missing BUSCOs
        429  |   Total BUSCO groups searched


### 4) LMP processing
There are two LMP libraries. To avoid repetition, we have taken some of these steps for one library only, but they should be performed with both.
 
a) Run FastQC to check read metrics for LMP as shown above.

b) Run the Python script to remove Nextera adapters from LMP reads and any PE contamination.  

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

The processed LMP files will be written to the 'nextclip' directory. These should be used in the subsequent scaffolding. The read counts before and after trimming are written to the log file, for the test dataset we get the following before and after values:

```
LIB6470 read count before trimming: 576252
LIB6471 read count before trimming: 576252

LIB6470 read count after trimming: 435974
LIB6471 read count after trimming: 443411
```

### 5) QC processed LMPs 
a) Use KAT comp to check for LMP representation issues by comparing LMP reads to PE reads to check for LMP representation issues 

```
kat comp -n -t 16 -m 27 -H10000000000 -I10000000000 -o lmp_vs_pe '/path/to/trimmed_lmp_R1_lib1.fastq /path/to/trimmed_lmp_R2_lib1.fastq' '/path/to/pe_R1.fastq /path/to/pe_R2.fastq'

kat plot spectra-mx -o lmp_vs_pe_spectra_mx.png -x 100 --intersection lmp_vs_pe-main.mx
```

<img src="images/lmp_vs_pe_spectra_mx.png"  width="450" height="400">

As they come from the same sample, there should be no content in the LMP data which is not present in the PE data, and vice versa. Hence, there should be no exclusive content.

b) Map the reads to a reference and generate an insert size histogram to check the insert size and the shape of the distribution. 

```
bwa index -p yeast ./contigs/a.lines.fasta
bwa mem -SP -t 8 yeast /path/to/trimmed_lmp_R1_lib1.fastq /path/to/trimmed_lmp_R2_lib1.fastq > lmp2contig.sam

bioawk -c'sam' '{if ($mapq>=60){if($tlen<0){print int($tlen/100)*100}else{print -int($tlen/100)*100}}}' lmp2ref.sam  | sort -n | uniq -c | awk '{print $2","$1}' > lmp_insert_sizes.txt

```

This is the expected distribution of the insert sizes of library 1: 

<img src="images/yeast_lmp.png"  width="500" height="400">

The distribution has a clear, pronounced peak so it is easy to see that the insert size is approximately 5000. There is no paired end contamination present, as this would cause another peak closer to the origin.  



c) Calculate the fragment coverage from trimmed read count and insert size

TODO

### 6) Scaffolding

s_scaff is an edited version of SOAPdenovo which is better suited to complex genomes, and more configurable.

a) Make a [SOAPdenovo config file] (http://soap.genomics.org.cn/soapdenovo.html) using both the PE and LMP reads to scaffold. 

```
[LIB]
avg_ins=320
q1=/path/to/pe_R1.fastq
q2=/path/to/pe_R2.fastq

[LIB]
avg_ins=7000
reverse_seq=1
q1=/path/to/trimmed_lmp_R1_lib1.fastq
q2=/path/to/trimmed_lmp_R2_lib1.fastq

[LIB]
avg_ins=15000
reverse_seq=1
q1=/path/to/trimmed_lmp_R1_lib2.fastq
q2=/path/to/trimmed_lmp_R2_lib2.fastq

```

The config file must be correctly configured, and there are lots of options to customize the configuration, details of which can be found in the [SOAP denovo documentation] (http://soap.genomics.org.cn/soapdenovo.html). It is advisable to familiarize yourself with these by varying them and observing the impact different parameters have on the final assembly.

We have kept our configuration file relatively simple, specifying only the paths to the data sets to be used for scaffolding, their type, and their insert size. The `reverse_seq` field indicates whether we have paired end (= 0) or long mate pair (= 1) read sets.
 
b) Run "prepare->map->scaff" pipeline.  


```
/path/to/s_prepare -t 8 -g yeast -K 71 -c ./contigs/a.lines.fasta 2>&1

PREFIX="yeast"
CONFIG_FILE="./soap.config"
NCPUS="32"

s_map -k 31 -s $CONFIG_FILE -p $NCPUS -g $PREFIX > $PREFIX.map.log 2>&1


s_scaff -p $NCPUS -g $PREFIX > $PREFIX.scaff.log 2>&1
```

The prepare script converts the fasta file to the format needed for SOAP mapping. For the majority of use cases, it is suitable to use K = 71. In s_map, K must be lower as reads may have been trimmed, and a lower value enables reads containing a small number of errors to be mapped to the contigs. 

Before proceeding from the map step to the scaffolding step, you should check that the mapping results are as expected, as if there are any problems at this stage, then the scaffolding step will not give good results. In particular, you should check that a reasonable proportion of reads have mapped to the contigs. The key part of the log for `s_map` for our example is:

```
Total reads         4623823
Reads in gaps       326317
Ratio               7.1%
Reads on contigs    3012707
Ratio               65.2%

```

The reads in gaps are reads which do not fall on a contig at all. As we have reasonable read coverage and have used a kat plot o check that information from the reads is not missing in the contigs, most of the reads should map at least partially to the contigs, which they do. The total number of reads is larger than the number of reads in gaps plus the number of reads on contigs because some reads map partially to a contig, and hang off the end. 

After the mapping has completed successfully, it's time to do the scaffolding. By default, `s_scaff` will scaffold in insert size order, from the smallest to the largest. To change this order, specify the `rank` field in the configuration file. Though `s_scaff` makes its own insert size estimates, it bases this ordering on the user specified insert sizes, so it is important that these are correct. You can check the following part of the log to make sure that the insert sizes calculated by `s_scaff` are similar to those specified in the config file:

 ```
 
 For insert size: 300
 Total PE links                      1141797
 Normal PE links on same contig      1131394
 Incorrect oriented PE links         1239
 PE links of too small insert size   7894
 PE links of too large insert size   0
 Correct PE links                    1089
 Accumulated connections             348
Use contigs longer than 300 to estimate insert size:
 PE links               1131017
 Average insert size    319
 SD                     138
 
 ```
 
If the `Aevrage insert size` field is missing, then the software has not been able to calculate it, which indiates that there is a serious issue with either the data or the configuration. As most contigs are significantly longer than the average insert size of the paired end library. most paired end reads map to the same contig. A PE link is a part of the graph which a read and its pair connect. Their distance apart on the graph must be consistent with the insert size. If two edges of the graph are linked by a read pair, then we have an 'Accumulated connection.'

The scaffolding begins after the reads have been loaded onto the graph.  

If this pipeline runs successfully, a number of output files will be created. The final scaffolds have the extension 'scafSeq.' 

c) SOAPdenovo converts gaps in contigs to Cs and Gs so we need to convert these back to Ns using the script included. The three input files are output by SOAP.

```
python SOAP_n_remapper.py <contigPosInScaff_file> <scafSeq_file> <contig_file> <output_file>
```

### 7) Scaffold validation
a) Check N50 and total content.  

* Total content: 11.78e6bp
* N50: 531425bp

The total content is similar to the expected genome size, so the assembly contains roughly the right amount of information. The N50 is reasonable for a genome of this size and complexity.

b) Use KAT comp to generate a spectra-cn to compare PE reads to scaffolds  

```
kat comp -t 16 -m 31 -H10000000000 -I10000000000 -o reads_vs_scaffolds '/path/to/pe_R1.fastq /path/to/pe_R2.fastq' /path/to/scaffolds/yeast.scafSeq
```

<img src="images/lmp_vs_pe-main.mx.spectra-cn.png"  width="450" height="400">

Again, there is no content from the reads missing in the assembly and no duplication of content, but there are a few erroneous kmers present.

c) Run QUAST and align sequences to BUSCOs. 

```
python /path/to/busco2/BUSCO.py -o busco_lmp -in ./yeast_ns_remapped.fasta -l ~/busco_data/eukaryota -m genome -f
```
   Count      |       Type
------------- | ------------------------------------
        417   |   Complete BUSCOs
        301   |   Complete and single-copy BUSCOs
        116   |   Complete and duplicated BUSCOs
        6     |   Fragmented BUSCOs
        6     |   Missing BUSCOs
        429   |   Total BUSCO groups searched

```
mkdir quast
python /path/to/quast/quast.py --extensive-mis-size 10000 -o ./quast -R ./yeast.scafSeq -t 8 -f ref/S288C_reference_sequence_R64-2-1_20150113.fsa
```
The `--extensive-mis-size` parameter sets a threshold for what is considered to be a local misassembly. By specifying this to be larger than the default value, we exclude very small rearrangements from the misassembly count. As the N50 of the contigs was greater than the average size of a gene, scaffolding did not increase the number of BUSCOs present. 

Genome statistics	 | yeast.scafSeq
-------------------- |---------------
Genome fraction (%)			  |	93.625
Duplication ratio			  |	1.260
Largest alignment			  |	414518
Total aligned length		  |	13152876
NGA50							  |	104027
LGA50							  |	31
Misassemblies					  |
misassemblies					  |71
Misassembled contigs length  |	4766054
Mismatches					  |
mismatches per 100 kbp		  | 4.10
indels per 100 kbp			  |1.78
N's per 100 kbp				  |	8291.35
Statistics without reference |	
contigs						  | 	279
Largest contig				  |	609133
Total length					  |	14473629
Total length (>= 1000 bp)	  | 14327082
Total length (>= 10000 bp)	  | 14176994
Total length (>= 50000 bp)	  | 11742048
Predicted genes	            |
predicted genes (unique)    |	7354

We can see that the scaffolder has successfully put together a large number of contigs without significantly increasing the number of misassemblies, which indicates that the scaffolds have been constructed correctly. 
