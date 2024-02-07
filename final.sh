##for Quality Control
fastqc 10610_1.fastq.gz
firefox 10610_1_fastqc.html

##for trimming the fastq files
java -jar /usr/src/Trimmomatic-0.36/trimmomatic-0.36.jar PE  \
-threads 4 -trimlog 10610.log \
10610_1.fastq.gz \
10610_2.fastq.gz \
10610_trimmed_R1_paired.fastq.gz \
10610_trimmed_R1_unpaired.fastq.gz \
10610_trimmed_R2_paired.fastq.gz \
10610_trimmed_R2_unpaired.fastq.gz \
ILLUMINACLIP:/usr/src/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36
##if error unable to detect quality score use either -phred33 or -phred64

##index the ref genome (chr1)
bwa index -a bwtsw shr.fasta

##map the reads to the ref genome using bwa mem
bwa mem -t 8 -R '@RG\tID:rg1\tSM:10610\tPL:illumina\tLB:lib1\t:PU:H522HDMXX:1:CCATCCTC+ATTACTCG' shr.fasta 10610_trimmed_R1_paired.fastq.gz 10610_trimmed_R2_paired.fastq.gz > 10610.sam

##Clean up and turning the SAM file to BAM
samtools fixmate -O bam 10610.sam 10610.bam

##to find the mapped reads
samtools view -F 4 10610.bam > mapped.sam

##to check the number of reads who are mapped fully
 cat mapped.sam | awk '{print $6}' | grep -c -P "^\d*M\b"

##Should index the genome assembly again(chr1)
samtools faidx shr.fasta

##create sequence directory
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar CreateSequenceDictionary \
      -R shr.fasta


##add read groups
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar AddOrReplaceReadGroups \
      -I 10610.bam \
      -O 10610_RG.bam \
      -RGID rg1 \
      -RGLB lib1 \
      -RGPL illumina \
      -RGPU H522HDMXX:1:CCATCCTC+ATTACTCG \
      -RGSM 10610

##sort by coordinate using picard
java -Djava.io.tmpdir=tmp -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SortSam \
     -I 10610.bam \
     -O 10610_sorted.bam \
     --CREATE_INDEX true \
     -SO coordinate

##mark dublicates using picard
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar MarkDuplicates \
      -I 10610_sorted.bam \
      -O 10610_sorted_marked.bam \
      --CREATE_INDEX true \
      -M 10610_sorted_marked_metrics.txt


##base recalibration using gatk
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar BaseRecalibrator \
   -I 10610_sorted_marked.bam \
   -R shr.fasta \
   --known-sites dbSNP.vcf.gz \
   -O 10610_sorted_marked_recal.table

##applying base recalibration
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar ApplyBQSR \
   -R shr.fasta \
   -I 10610_sorted_marked.bam \
   --bqsr-recal-file 10610_sorted_marked_recal.table \
   -O 10610_sorted_marked_bqsr.bam

##haplotype caller to get GVCF file
java -Xmx4g -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller  \
   -R shr.fasta \
   -D dbSNP.vcf.gz \
   -I 10610_sorted_marked_bqsr.bam \
   -O 10610_aln_gvcf.g.vcf.gz \
   -ERC GVCF



##genotypeVcf to get vcf from gvcf
java -Xmx4g -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GenotypeGVCFs \
  -R shr.fasta \
  -D dbSNP.vcf.gz \
  -V 10610_aln_gvcf.g.vcf.gz \
  -O 10610_aln_raw.vcf.gz

 ##selecting SNPs from raw vcf
 java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
      -R shr.fasta \
      -V 10610_aln_raw.vcf \
      --select-type-to-include SNP \
      -O 10610_aln_rawSNP.vcf.gz

##to get those of the indels
##java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar SelectVariants \
##     -R shr.fasta \
##     -V 10610_aln_raw.vcf \
##     --select-type-to-include INDEL \
##     -O 10610_aln_rawIndel.vcf.gz

##Filter SNPs
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
   -R shr.fasta \
   -V 10610_aln_rawSNP.vcf.gz \
   -O 10610_aln_filterSNP.vcf.gz \
   --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0" \
   --filter-name "SNPbasic"

##Filter SNPs genotype
java -jar gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar VariantFiltration \
-R shr.fasta \
-V 10610_aln_filterSNP.vcf.gz \
-G-filter "GQ < 20" \
-G-filter-name "lowSNPGQ" \
-O 10610_aln_filterGQSNP.vcf.gz
