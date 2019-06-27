# TILLinG-mutants
Bioinformatics pipeline for the analysis of TILLinG mutants in conjunction with the new optimised reference sequence for tetraploid wheat
<br/>Pipeline Commands:
<br/><b>Mutation Calling Pipeline</b></br>
A) Generate pileup file for 100 samples:</br>
>samtools mpileup  -B -A -d 99999 -q 10  -f parent_derived_ref.fa  -t DP -b bam_files_100_samples.txt > samples.pileup</br>
Where -q 10 consider only reads with MAPQ of 10 or more. bam_files_100_samples.txt contains name of sorted bam files for parental and 99 mutant samples.</br>

B)Mutation Calling and generate Mutation file:
>java -jar MFbio.jar showform=no task=mutation2  srcdir=samples.pileup.gz destdir=/dierectoy  p1=100  p2=3  p3=0.1
Where srcdir is pileup file generated in step A; destdir is directory where output file unique_mutation.mut will be generated ; p1 is number
of samples in pileup file; p2 is minimum number of reads for a non-reference allele (3 recommended); p3 is minimum ratio for minor allele 
to call a position heterozygous. For example, if p3=0.1 and allele 1 has 95% and allele 2 has 5% of reads then position is called homozygous,
whereas if allele 1 has 85% and allele 2 has 15% of reads then position is called heterozygous.

<b>Discovery of highly covered genomic regions Pipeline:</b></br>
A)Parental sample were mapped to Kronos genomic reference and a pileup file was generated.</br>  
samtools mpileup  -B -A -d 99999 -q 10 -f kronos_EI_v1.fasta  -t DP -b parent_bam_file.txt > parent_kronos.pileup</br>
Where -q 10 consider only reads with MAPQ of 10 or more. parent_bam_file.txt contains name of sorted bam file for parental sample.</br>
B)Genomic regions that were significantly covered were extracted.</br>
>java -jar  MFbio.jar showform=no task=highcoverage srcdir=parent_kronos.pileup.gz file1=kronos_EI_v1.fasta.fai  destdir=kronos_high_coverage_17.bed p1=17 p2=500 p3=300</br>
Where srcdir is generated pileup file in step A; file1 is index of Kronos fasta file; destdir refers to generated output bed file; p1 minimum depth or coverage; p2 tail length; p3 maximum distance to merge 2 coninuous regions.




