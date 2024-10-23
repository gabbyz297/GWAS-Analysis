_See the code this pipeline is derived from [here](https://github.com/CassinSackett/SNP_capture)_

## Run FastQC on all files to check quality of reads :white_check_mark:
`FastQC` is a program used to check the sequence quality of your reads. This information can be used to determine if a sample needs to be resequenced or if bad sequence quality base pairs can be removed.

`for i in /path/to/files/*.fq do /path/to/FastQC/fastqc $i; done`

_To see how to run this code in parallel click [here](https://github.com/gabbyz297/GWAS-Analysis/blob/master/fastqc.sh)_

## Trim reads using Trimmomatic :scissors:
`Trimmomatic` trims reads based on the average quality score within a sliding window. The sliding window parameters are set by using the first number as your desired size of the window and the second number as the minimum quality score to keep the base pairs within that window. You can specify a number of leading or trailing base pairs to trimmed automatically based on `FASTQC` results or you can use the default settings below. 

_1P= paired read 1, 1U= unpaired read 1, 2P= paired read 2, 2U= unpaired read 2_

`module load jdk/1.8.0_262/intel-19.0.5`

`java -jar /path/to/Trimmomatic-0.39/trimmomatic-0.39.jar PE /path/to/forward/file.fq /path/to/reverse/file.fq /path/to/forward/file/file_1P.fq /path/to/forward/file/file_1U.fq /path/to/reverse/file/file_2P.fq /path/to/reverse/file/file_2U.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

## Index Reference Genome using BWA Index; Create .fai and .dict files with Samtools and GATK :open_book:
Indexing the genome is like a book index, used to make alignment easier later 

### BWA Index

_`-a` algorithm to construct index from, `bwtsw` algorithm implemented in BWT-SW which is for large genomes, just need `-a` for small genomes_

`/path/to/bwa index -a bwtsw /path/to/reference/reference.fa`

### Samtools

`module load samtools/1.10/intel-19.0.5
samtools faidx /path/to/reference/reference.fa`

### GATK

_`-Xmx2G` sets the initial and maximum heap size available to improve performance_

`/path/to/gatk-4.1.2.0/gatk --java-options "-Xmx2G" CreateSequenceDictionary -R /path/to/reference/reference.fa`


## Align reads to reference using BWA mem :left_right_arrow:
Now that the reference genome is indexed we can align our reads to it using `BWA mem`

_`bwa mem` aligns 70bp-1Mbp sequences, `-v` verbose level (a value 0 for disabling all the output to stderr; 1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging), `-M` mark shorter split hits as secondary (for Picard compatibility), `-P` paired-end mode (only necessary if you have your forward and reverse reads in one interleaved file), `-a` output all found alignments, `-t` number of threads_

**outputs end in 1P.sam but contains both read 1 and 2**

`module load bwa/0.7.17/intel-19.0.5`

`for i in /path/to/files/*1P.fq; do bwa mem -a -M -t 10 -v 3 \ 
/path/to/reference/reference.fa \ 
$i ${i/1P/2P} > ${i%P.fq}.sam 2> ${i%P.fq}.mem.log; done`

## Create BAM files using Samtools View 👓
Next we will convert our SAM files into binary BAM files to speed up downstream analyses and save storage space.

_`-q` skip alignments with mapping quality less than specified number, -bt output in BAM format_

`module load samtools/1.10/intel-19.0.5`

`for i in /path/to/files/*.sam; do samtools view -q 20 -bt /path/to/reference/.fa -o ${i%sam}bam $i; done`

## Check Alignment Statistics with Samtools :white_check_mark:
BWA mem doesn't ouput it's own alignment stats so I've added this step as a checkpoint before filtering and genotyping.

`module load samtools/1.10/intel-19.0.5`

`for i in /path/to/files/*.bam; do samtools flagstat -o ${i%bam}bam_stat $i; done`

## Assign all reads to sorted read-groups using Picard tools :books:
Next we will assign our reads to read-groups based on which lane of the sequencer they were run on. This is done to minimize the impact of sequencer effects on results.

_`-Dpicard.useLegacyParser=false` for GATK 4.0 to run Picard tools from GATK, `-Xmx2G` sets the initial and maximum heap size available to improve performance,`AddOrReplaceReadGroups` assigns reads to read group, `-I` input file, `-O` output file, `-MAX_RECORDS_IN_RAM` specifies the number of records stored in RAM before storing in disk, `TMP_DIR $PWD/tmp` creates a temorary directory called tmp, `-SO` sort order to ouput in, `-ID` read-group id, `-LB` read group library, `-PL` read group platform, `-PU` read group platform unit, `-SM` read group sample name_

`cd $PBS_O_WORKDIR
mkdir -p tmp`

`module load jdk/1.8.0_262/intel-19.0.5`

`for i in /path/to/files/*.bam;
do
        java -Dpicard.useLegacyParser=false -Xmx2g -jar /path/to/picard/picard.jar \
        AddOrReplaceReadGroups -I $i -O ${i%.bam}.tag.bam -MAX_RECORDS_IN_RAM 1000000 -TMP_DIR $PWD/tmp \
        -SO coordinate -ID ${i%.bam} -LB 1 -PL illumina -PU 1 -SM ${i%.bam};
done`

## Mark and remove PCR duplicates using Picard tools :wastebasket:

_`-cd` change directory, `mkdir -p` make directory and if neccessary, all parent directories, `-Xmx2G` sets the initial and maximum heap size available to improve performance, `MarkDuplicates` identifies duplicate reads, `MAX_FILE_HANDLES_FOR_READ_ENDS_MAP` maximum number of file handles to keep open when spilling read ends to disk (default is 8000), `MAX_RECORDS_IN_RAM` when writing files that need to be sorted this will specify the number of records stored in RAM before spilling to disk, `TMP_DIR $PWD/tmp` creates a temorary directory called `tmp`, `METRICS_FILE File` to write duplication metrics to, `ASSUME_SORTED` If true, assume that the input file is coordinate sorted_

`module load  jdk/1.8.0_262/intel-19.0.5
cd $PBS_O_WORKDIR
mkdir -p temp`

`for i in /path/to/files/*.tag.bam; 
do
        java -Xmx4g -jar /path/to/picard/picard.jar \
        MarkDuplicates -INPUT $i -OUTPUT ${i%.tag.bam}.rmdup.bam -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 6000 -MAX_RECORDS_IN_RAM 1000000 -TMP_DIR $PWD/temp \
        -METRICS_FILE ${i%.tag.bam}.rmdup.metrics -ASSUME_SORTED true;
done`

## Index filtered files using Samtools :open_book:
Next filtered files will be indexed for realignment with `samtools`

`module load samtools/1.10/intel-19.0.5`

`for i in /path/to/files/*.rmdup.bam; do samtools index $i; done`

## Create GVCF files using GATK :file_folder:
Next we will convert our filtered BAM files to GVCF files which will create VCF formatted files without genotypes

_`-Xmx2G` sets the initial and maximum heap size available to improve performance, HaplotypeCaller `--ERC GVCF` creates GVCF, `--min-base-quality-score` minimum quality required to consider a base for calling_

`path/to/gatk-4.1.2.0/gatk --java-options "-Xmx60g" HaplotypeCaller --ERC GVCF -R /path/to/reference/.fa  -I /path/to/file/file.rmdup.bam -O /path/to/file/file.g.vcf --min-base-quality-score 20`

## Combine GVCF files using GATK Combine GVCF :heavy_plus_sign:
We can combine samples into a single GVCF file for genotyping.

`path/to/gatk-4.1.2.0/gatk CombineGVCFs -R /path/to/reference/.fa  -V /path/to/file1/file1.g.vcf -V /path/to/file2/file2.g.vcf -O /path/to/file/combined.g.vcf` 

## Genotype GVCF file :dna:
Now that we have a single GVCF file we can genotype it.

`path/to/gatk-4.1.2.0/gatk --java-options "-Xmx60g" GenotypeGVCFs -R /path/to/reference/.fa  -V /path/to/file/combined.g.vcf -O /path/to/file/combined.vcf`

## Select SNP variants using GATK :magnet:
Now that we have a genotyped VCF file we can select variant types that we're interested in. 
 _`--select-type-to-include` options: INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION_
 
`/path/to/gatk-4.1.2.0/gatk SelectVariants --variant /path/to/file/combined.vcf -R /path/to/reference/.fa --output /path/to/file/file2.vcf --select-type-to-include SNP --exclude-non-variants true --set-filtered-gt-to-nocall true` 

## Quality filter variants and flag variants that don't meet criteria using GATK :wastebasket:
Next we can filter our genotyped VCF file for analyses.

_`ReadPosRankSum` compares whether positions of reference and alternate alleles are different within reads, `MQRankSum` compares mapping qualities of reads supporting the reference allele and alternate allele, `FS` Phred-scaled probability that there is strand bias at the site. `FS` value will be close to 0 when there is little to no strand bias at the site, `QD` is the variant confidence divided by the unfiltered depth of samples. This normalizes the variant quality in order to avoid inflation caused when there is deep coverage, `DP` is genotype depth of coverage_

`/path/to/gatk-4.1.2.0/gatk VariantFiltration --variant /path/to/file/file2.vcf --output /path/to/file/file3.vcf -R /path/to/reference/.fa \
--filter-name "ReadPosRankSum_filter" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "MQRankSum_filter" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "FS_filter" \
--filter-expression "FS > 60.0" \
--filter-name "QD_filter" \
--filter-expression "QD < 2.0" \
--genotype-filter-name "DP8filter" \
--genotype-filter-expression "DP < 8" 2>/dev/null` 

## Remove flagged variants using GATK :wastebasket:
This will remove variants that did not meet the filters specified above.

`module load jdk/1.8.0_262/intel-19.0.5`

`/path/to/gatk-4.1.2.0/gatk SelectVariants -R /path/to/reference/.fa --variant /path/to/file/file3.vcf --output /path/to/file/file4.vcf --set-filtered-gt-to-nocall true`

## View general stats of VCF file using VCFtools 👓

`Vcftools` can output a variety of stats including number of SNPS, number of heterozygous sites, number of homozygous sites etc.

_`PERL5LIB` allows you to use perl through `vcftools`_

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcf-stats /path/to/file/.vcf  > output.txt`

## Remove poor sequences with VCFtools 🗑️
The first part of this code outputs a file which has 5 columns: `N_MISS`= number of sites the individual does not have data for, `F_MISS`= frequency of missing data for the individual

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --missing-indv --out /path/to/file/file`

This ouput can be used to create a text file of individuals to remove based on missing data. This file will be used in the code below. This step can be skipped if no individuals have a high percentage of missing data.

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --remove file.txt --out /path/to/file.vcf`

## Create chromosome mapping file using bcftools view for VCF to PLINK conversion 🧬
List all the chromosomes or contigs in your `VCF` file, one chromosome or scaffold per line for `PLINK`. 

`/path/to/bin/bcftools view -H file.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > /path/to/output/chrom-map.txt`

## Convert VCF to PLINK (.ped + .map) file for analyses (vcftools) :arrows_counterclockwise:
This will create files compatible with PLINK 1.9

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file/file.vcf --chrom-map  /path/to/file/file.txt --out /path/to/file/output  --plink`

## Convert PLINK (.ped + .map) file to PLINK2 (.pgen + .psam + .pvar) (plink) :arrows_counterclockwise:
This will create files compatible with PLINK 2.0

`/path/to/plink2 --pedmap /path/to/prefix_of_ped_and_map/file --allow-extra-chr --out output_prefix`

## Create phenotype file for PLINK association test 📁
We need to create a phenotype files so `PLINK` knows how to group our samples for our GWAS. The text file format uses FULL NAME from `VCF` file. My `VCF` file named all my individuals the full path (pretty sure that's the default) to the file so, double check what your individual names are in the `VCF` file you will be using for your association test.

_The phenotype files has a 3 column format with Family ID (`FID`), Individual ID (`IID`) and phenotype [1= unaffected (or can be control), 2= affected (or can be case), 0= missing data], tab delimited columns. If there is no Family ID you can copy the Individual ID into that column or remove it completely (I haven't tried that though)._

`path/to/file    path/to/file    1,2 or 0`

## Run PCA in PLINK 2.0 🏃‍♀️
We can run a PCA to make sure that the results of our GWAS are not due to any population structuring

_`--allow-extra-chr` allow extra chromosomes [permit unrecognized chromosome codes, treating their variants as unplaced]_

`/path/to/plink2 --pfile /path/to/prefix/file --allow-extra-chr --pheno /path/to/phenotype/file.txt --pca --out output`

## Calculate Allele Frequencies in PLINK 2.0 (for small sample sizes) 🧬
  This creates a table with allele frequencies output as a `.afreq` file, This can be used to account for differences in minor allele frequencies when running an associaation test 
 
`/path/to/plink2 --pfile /path/to/prefix/file --allow-extra-chr --pheno /path/to/phenotype/file.txt --freq --out /path/to/prefix/file`

## Run association test in PLINK 2.0 🏃‍♀️
We can include our top PCs in our association test as a covariate using `--covar` and specify the number of PCs with `--parameters` to account for population structure. 

`/path/to/plink2 --pfile /path/to/prefix/file --allow-extra-chr --pheno /path/to/phenotype/file.txt --read-freq /path/to/file.afreq --glm --covar /path/to/file.eigvec  --parameters 1-4 --out /path/to/prefix/file`


## Additional Analyses 📊

### Calculate Tajima's D with Vcftools 

  _`--TajimaD` specify sliding window size for Tajima's D to be calculated_
  
`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file.vcf --out /path/to/file --TajimaD 10000`



### Calculate Fst with Vcftools
Fst can be calculated between two populations with `vcftools` by creating two text files with the names of the samples from each group as they appear in the vcf file.

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file.vcf --weir-fst-pop population_1.txt --weir-fst-pop population_2.txt --out /path/to/file`

### Calculate Linkage Disequilibrium with Vcftools

`export PERL5LIB=/path/to/vcftools/src/perl
/path/to/vcftools/bin/vcftools --vcf /path/to/file.vcf --out /path/to/file --geno-r2 --ld-window-bp 50000`










