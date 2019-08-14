
# TILLinG-mutants
Bioinformatics pipeline for the analysis of TILLinG mutants in conjunction with the new optimised reference sequence for tetraploid wheat.

## Background
The pipeline was developed as part of a joint project between Australian Centre for Plant Functional Genomics (ACPFG), University of Saskatchewan and Università degli Studi della Tuscia for the project “A Molecular Diversity Drive for Precision-Engineered Wheat”.

Seeds of an advanced durum breeding line, UAD0951096_F2:5 (which has given rise to the commercially grown cultivar DBA-Aurora), were mutagenized with 0.7% EMS and 500 mutant plants were grown through one further generation to seed. In the next generation (M2) DNA from 99 of these plants and the unmutagenized control (UAD0951096_F2:5) were sequenced (Illumina, 2x150 bp) following reduction by exome capture using the 106.9 Mb [Roche NimbleGen SeqCap](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/shareddesigns.html) Wheat Exome Design. 

To improve variant calling, a novel consolidated Durum Exome Capture Reference (DECaR) was built: reads from the unmutagenized control sample were mapped into two publically available durum wheat whole genome assemblies of [Svevo](https://www.interomics.eu/durum-wheat-genome) and [Kronos EI v1](https://opendata.earlham.ac.uk/opendata/data/Triticum_turgidum/EI/v1.1/), and then combined with assembled unmapped control reads. The resulting reference was 484.5 Mbp.

Reads of the exome-captured samples were aligned to DECaR. Any variation from DECaR was considered a mutation if it was present in only one mutant sample and non-polymorphic between the control and DECaR. We also demanded a mutation to be covered by at least three reads, and had sufficient coverage for the control allele in at least 50 other mutant samples. DECaR also allowed us to predict zygosity of the mutations discovered.

All mutations were deposited into a database that is publicly available at http://duwtill.acpfg.com.au.

For further details and to cite:
>Mario Fruzangohar, Elena Kalashyan, Priyanka Kalambettu, Jennifer Ens, Krysta Wiebe, Curtis J. Pozniak, Penny J. Tricker, Ute Baumann "Novel Informatics Tools to Support Functional Annotation of the durum wheat genome." Frontiers in Plant Science (In review)

## Pipeline Commands
To be able to run the pipeline, make sure Java 1.8 or later is installed on your machine. Then download `WHEATbio.jar` from this repository and simply refer to the `WHEATbio.jar` file in all the following commands.
#### Discovery of highly covered genomic regions to construct DECaR
A. Map the control sample to the genomic reference, for example Kronos, and generate a pileup file.
```bash  
samtools mpileup -B -A -d 99999 -q 10 -f kronos_EI_v1.fa -t DP -b control_bam_file.txt > control_kronos.pileup
```
Where `-q 10` considers only reads with MAPQ of 10 or more and `parent_bam_file.txt` contains the name of the sorted bam file for the control sample.
B. Extract genomic regions that were significantly covered.
```bash
java -jar WHEATbio.jar showform=no task=highcoverage srcdir=control_kronos.pileup.gz file1=kronos_EI_v1.fa.fai destdir=kronos_high_coverage_17.bed p1=17 p2=500 p3=300
```
Where `srcdir` is the generated pileup file in step A; `file1` is the index of the Kronos fasta file; `destdir` refers to the generated output bed file; `p1` is the minimum depth or coverage; `p2` is the tail length and `p3` is the maximum distance to merge two continuous regions.
#### Mutation Calling
A. Generate the pileup file for all 100 samples.
```bash  
samtools mpileup -B -A -d 99999 -q 10 -f DECaR.fa -t DP -b bam_files_100_samples.txt > samples.pileup
```
Where `-q 10` considers only reads with MAPQ of 10 or more and `bam_files_100_samples.txt` contains the names of sorted bam files for both the control and 99 mutant samples.

B. Call the mutations.
```bash  
java -jar WHEATbio.jar showform=no task=mutations srcdir=samples.pileup.gz destdir=/directoy p1=100 p2=3 p3=0.1
```
Where `srcdir` is the pileup file generated in step A; `destdir` is the directory where the output file `unique_mutation.mut` will be generated; `p1` is the number of samples in the pileup file; `p2` is the minimum number of reads for a non-reference allele (3 is recommended); `p3` is the minimum proportion for a minor allele to call a position heterozygous. For example, if `p3=0.1` and allele 1 has 95% and allele 2 has 5% of reads then the position is called homozygous, whereas if allele 1 has 85% and allele 2 has 15% of reads then the position is called heterozygous.
#### Extraction of flanking sequences for mutation positions in the control sample
A. Map the control sample to DECaR and generate the pileup file.
```bash
samtools mpileup -B -A -d 99999 -q 10 -f DECaR.fa -t DP -b control_bam_file.txt > control_DECaR.pileup
```
B. Generate the variant file for the control file.
```bash
java -jar WHEATbio.jar showform=no task=variant srcdir=control_DECaR.pileup.gz destdir=control.var p1=1 p2=2 p3=0.1
```
Where `srcdir` is the pileup file generated at step A; `destdir` is the name of the output variant file; `p1` is the number of samples in the pileup file; `p2` and `p3` are the same as in the mutation calling step above (`task=mutations`).

C. Generate the flanking sequences of mutation positions.
```bash
java -jar WHEATbio.jar showform=no task=buildsampleseq srcdir=/directoy/unique_mutation.mut p1=50 p2=200 1=DECaR.fa 2=DECaR.fa.fai file1=control.var
```
Where `srcdir` is the mutation file generated by the mutation calling pipeline; `p1` is the minimum length of the flanking sequence; `p2` is the maximum length of the flanking sequence; `1` is the DECaR fasta file; `2` is the index of the DECaR fasta file; `file1` is the variant file of the control sample generated at step B.

