# mansonella

<!-- MarkdownTOC levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
	- [Software](#software)
	- [Directories](#directories)
	- [Create directories](#create-directories)
- [Assemble mitochondrial genome of Mansonella using 3 different strains](#assemble-mitochondrial-genome-of-mansonella-using-3-different-strains)
	- [Download Mansonella ozzardi mitochondra reference](#download-mansonella-ozzardi-mitochondra-reference)
	- [Align Illumina reads from each strain to Mansonella ozzadi mitochondria using BWA MEM](#align-illumina-reads-from-each-strain-to-mansonella-ozzadi-mitochondria-using-bwa-mem)
	- [Extract reads that mapped to the Mansonella ozzadi mitochondria reference](#extract-reads-that-mapped-to-the-mansonella-ozzadi-mitochondria-reference)
	- [Align MinION reads from each strain to Mansonella ozzadi mitochondria](#align-minion-reads-from-each-strain-to-mansonella-ozzadi-mitochondria)
	- [Assemble Mansonella mitochondrial genome with Unicycler](#assemble-mansonella-mitochondrial-genome-with-unicycler)
		- [Assemble T7_2 and T-8 using long read FASTQ and subset short read FASTQs](#assemble-t7_2-and-t-8-using-long-read-fastq-and-subset-short-read-fastqs)
		- [Assemble T6 using subset short read FASTQs](#assemble-t6-using-subset-short-read-fastqs)
	- [Create contig files for NOVOPlasty assembly](#create-contig-files-for-novoplasty-assembly)
	- [Assemble Mansonella mitochondrial genome using short read FASTQs with NOVOplasty](#assemble-mansonella-mitochondrial-genome-using-short-read-fastqs-with-novoplasty)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
PERL_BIN_DIR=/usr/local/packages/perl-5.24.0/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

BWA_BIN_DIR=/usr/local/packages/bwa-0.7.17/bin
NOVOPLASTY_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/NOVOPlasty_v3.8.1
PILON_BIN_DIR=/usr/local/packages/pilon-1.22
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
UNICYCLER_BIN_DIR=/usr/local/packages/python-3.5.2/bin




JULIA_BIN_DIR=/usr/local/bin
JULIA_LIB_DIR=/home/mattchung/.julia
PYTHON_BIN_DIR=/usr/local/packages/python-3.5.2/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

BBTOOLS_BIN_DIR=/usr/local/packages/bbtools-38.47
DOS2UNIX_BIN_DIR=/usr/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin

EXPRESS_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/express-1.5.1
FADU_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/FADU
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
KALLISTO_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/kallisto_v0.46.1
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
SUBREAD_BIN_DIR=/usr/local/packages/subread-1.6.4/bin
```

## Directories

```{bash, eval = F}
REFERENCES_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella/references
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella
OUTPUT_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella/output
```

## Create directories

```{bash, eval = F}
mkdir -p "$WORKING_DIR"/assemblies/unicycler/assembly1/T6
mkdir -p "$WORKING_DIR"/assemblies/unicycler/assembly1/T7_2
mkdir -p "$WORKING_DIR"/assemblies/unicycler/assembly1/T8

mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T6
mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T7_2
mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T8
```

# Assemble mitochondrial genome of Mansonella using 3 different strains

T6 = 5 larvae  
T8 = 3 larvae/adults  
T7.2 = 1 adult  

## Download Mansonella ozzardi mitochondra reference
```{bash, eval = F}
wget 
```

## Align Illumina reads from each strain to Mansonella ozzadi mitochondria using BWA MEM

##### Input Sets:
```{bash, eval = F}
REF_FNA=/local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
SEED_LENGTH=23
THREADS=16

## T6
OUTPUT_PREFIX=T6
FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz

## T7_2
OUTPUT_PREFIX=T7_2
FASTQ1=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R2_trimmed.fastq.gz

## T8
OUTPUT_PREFIX=T8
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
```

##### Commands:
```{bash, eval = F}
"$BWA_BIN_DIR"/bwa index "$REF_FNA"
echo -e ""$BWA_BIN_DIR"/bwa mem -t "$THREADS" -k "$SEED_LENGTH" "$REF_FNA" "$FASTQ1" "$FASTQ2" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".map1.bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N bwa -wd "$WORKING_DIR"/assemblies/
```

## Extract reads that mapped to the Mansonella ozzadi mitochondria reference

##### Inputs:
```{bash, eval = F}
## T6
OUTPUT_PREFIX=T6
FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz
BAM=/local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/T6.map1.bam

## T7_2
OUTPUT_PREFIX=T7_2
FASTQ1=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R2_trimmed.fastq.gz
BAM=/local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/T7_2.map1.bam

## T8
OUTPUT_PREFIX=T8
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
BAM=/local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/T8.map1.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F 4 "$BAM" | awk '{print $1}' | sort -n | uniq > "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".map1.mappedreads.list
"$SEQTK_BIN_DIR"/seqtk subseq "$FASTQ1" "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".map1.mappedreads.list | gzip > "$WORKING_DIR"//assemblies/"$(basename "$FASTQ1" | sed "s/[.]fastq.gz/.subset1.fastq.gz/g")"
"$SEQTK_BIN_DIR"/seqtk subseq "$FASTQ2" "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".map1.mappedreads.list | gzip > "$WORKING_DIR"/assemblies/"$(basename "$FASTQ2" | sed "s/[.]fastq.gz/.subset1.fastq.gz/g")"
```

## Align MinION reads from each strain to Mansonella ozzadi mitochondria
```{bash, eval = F}

```

## Assemble Mansonella mitochondrial genome with Unicycler

### Assemble T7_2 and T-8 using long read FASTQ and subset short read FASTQs

T6 was never sequenced on the MinION

##### Inputs:
```{bash, eval = F}
THREADS=16


## T7_2
OUTPUT_PREFIX=T7_2
LONG_FASTQ=
SHORT_FASTQ1=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R1_trimmed.fastq.gz
SHORT_FASTQ2=/local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R2_trimmed.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/assemblies/unicycler/assembly1/T7_2

## T8
OUTPUT_PREFIX=T8
LONG_FASTQ=/local/projects/RDMIN/SEQUENCE/20180806/fastq/emans_T8_mda_branched.fastq
SHORT_FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
SHORT_FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/assemblies/unicycler/assembly1/T8
```

##### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$UNICYCLER_BIN_DIR"/unicycler --mode normal --long "$LONG_FASTQ" --short1 "$SHORT_FASTQ1" --short2 "$SHORT_FASTQ2" -o "$OUTPUT_DIR" --pilon_path "$PILON_BIN_DIR"/pilon-1.22.jar -t "$THREADS"" | qsub -P jdhotopp-lab -N unicycler -wd "$OUTPUT_DIR" -q threaded.q -pe thread "$THREADS" -l mem_free=5G
```

### Assemble T6 using subset short read FASTQs
##### Inputs:
```{bash, eval = F}
THREADS=16

## T6
OUTPUT_PREFIX=T6
SHORT_FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
SHORT_FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/assemblies/unicycler/assembly1/T6
```

##### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$UNICYCLER_BIN_DIR"/unicycler --mode normal --short1 "$SHORT_FASTQ1" --short2 "$SHORT_FASTQ2" -o "$OUTPUT_DIR" --pilon_path "$PILON_BIN_DIR"/pilon-1.22.jar -t "$THREADS"" | qsub -P jdhotopp-lab -N unicycler -wd "$OUTPUT_DIR" -q threaded.q -pe thread "$THREADS" -l mem_free=5G
```

## Create contig files for NOVOPlasty assembly

T6: "$WORKING_DIR"/assemblies/novoplasty/assembly1/T6/config.txt
```{bash, eval = F}
Project:
-----------------------
Project name          = T6
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Extend seed directly  = no
Reference sequence    = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Variance detection    =
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = /local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
Reverse reads         = /local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3
Use Quality Scores    = no




Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Type                 = (chloro/mito/mito_plant) "chloro" for chloroplast assembly, "mito" for mitochondrial assembly and
                       "mito_plant" for mitochondrial assembly in plants.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 33).
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23.
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                  
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Extend seed directly = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed
                       originates from the same sample and there are no possible mismatches (yes/no)
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast.
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file
                       with all the variances compared to the give reference (yes/no)
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina/ion - The performance on Ion Torrent data is significantly lower
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)

Heteroplasmy:
-----------------------
MAF                  = (0.007-0.49) Minor Allele Frequency: If you want to detect heteroplasmy, first assemble the genome without this option. Then give the resulting
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option
                       (0.01 will detect heteroplasmy of >1%)
HP exclude list      = Option not yet available
PCR-free             = (yes/no) If you have a PCR-free library write yes

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.9).
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.3).
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the
                       300 bp reads of Illumina (yes/no)

```

T7_2: "$WORKING_DIR"/assemblies/novoplasty/assembly1/T7_2/config.txt
```{bash, eval = F}
Project:
-----------------------
Project name          = T7_2
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Extend seed directly  = no
Reference sequence    = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Variance detection    =
Chloroplast sequence  = 
Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = /local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R1_trimmed.fastq.gz
Reverse reads         = /local/projects/EMANS/T7_2/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106387_S2_L001_R2_trimmed.fastq.gz

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3
Use Quality Scores    = no




Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Type                 = (chloro/mito/mito_plant) "chloro" for chloroplast assembly, "mito" for mitochondrial assembly and
                       "mito_plant" for mitochondrial assembly in plants.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 33).
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23.
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                  
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Extend seed directly = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed
                       originates from the same sample and there are no possible mismatches (yes/no)
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast.
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file
                       with all the variances compared to the give reference (yes/no)
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina/ion - The performance on Ion Torrent data is significantly lower
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)

Heteroplasmy:
-----------------------
MAF                  = (0.007-0.49) Minor Allele Frequency: If you want to detect heteroplasmy, first assemble the genome without this option. Then give the resulting
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option
                       (0.01 will detect heteroplasmy of >1%)
HP exclude list      = Option not yet available
PCR-free             = (yes/no) If you have a PCR-free library write yes

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.9).
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.3).
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the
                       300 bp reads of Illumina (yes/no)

```

T8: "$WORKING_DIR"/assemblies/novoplasty/assembly1/T8/config.txt
```{bash, eval = F}
Project:
-----------------------
Project name          = T8
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Extend seed directly  = no
Reference sequence    = /local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
Variance detection    =
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = 151
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = /local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
Reverse reads         = /local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3
Use Quality Scores    = no




Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Type                 = (chloro/mito/mito_plant) "chloro" for chloroplast assembly, "mito" for mitochondrial assembly and
                       "mito_plant" for mitochondrial assembly in plants.
Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.
                       Default value for mito: 12000-20000 / Default value for chloro: 120000-200000
                       If the expected size is know, you can lower the range, this can be useful when there is a repetitive
                       region, what could lead to a premature circularization of the genome.
K-mer                = (integer) This is the length of the overlap between matching reads (Default: 33).
                       If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23.
                       For reads longer then 101 bp, this value can be increased, but this is not necessary.
Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited                  
                       memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB
                       (if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))
Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Extend seed directly = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed
                       originates from the same sample and there are no possible mismatches (yes/no)
Reference (optional) = If a reference is available, you can give here the path to the fasta file.
                       The assembly will still be de novo, but references of the same genus can be used as a guide to resolve
                       duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast.
                       References from different genus haven't beeen tested yet.
Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file
                       with all the variances compared to the give reference (yes/no)
Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).
                       You have to assemble the chloroplast before you assemble the mitochondria of plants!

Dataset 1:
-----------------------
Read Length          = The read length of your reads.
Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.
Platform             = illumina/ion - The performance on Ion Torrent data is significantly lower
Single/Paired        = For the moment only paired end reads are supported.
Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)
Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)
Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)

Heteroplasmy:
-----------------------
MAF                  = (0.007-0.49) Minor Allele Frequency: If you want to detect heteroplasmy, first assemble the genome without this option. Then give the resulting
                       sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option
                       (0.01 will detect heteroplasmy of >1%)
HP exclude list      = Option not yet available
PCR-free             = (yes/no) If you have a PCR-free library write yes

Optional:
-----------------------
Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)
Insert Range         = This variation on the insert size, could lower it when the coverage is very high or raise it when the
                       coverage is too low (Default: 1.9).
Insert Range strict  = Strict variation to resolve repetitive regions (Default: 1.3).
Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the
                       300 bp reads of Illumina (yes/no)

```

## Assemble Mansonella mitochondrial genome using short read FASTQs with NOVOplasty

##### Input Sets:
```{bash, eval = F}
## T6
CONFIG="$WORKING_DIR"/assemblies/novoplasty/assembly1/T6/config.txt
OUTPUT_DIR="$WORKING_DIR"/assemblies/novoplasty/assembly1/T6/

## T7_2
CONFIG="$WORKING_DIR"/assemblies/novoplasty/assembly1/T7_2/config.txt
OUTPUT_DIR="$WORKING_DIR"/assemblies/novoplasty/assembly1/T7_2/

## T8
CONFIG="$WORKING_DIR"/assemblies/novoplasty/assembly1/T8/config.txt
OUTPUT_DIR="$WORKING_DIR"/assemblies/novoplasty/assembly1/T8/
```

```{bash, eval = F}
echo -e ""$PERL_BIN_DIR"/perl "$NOVOPLASTY_BIN_DIR"/NOVOPlasty3.8.1.pl -c "$CONFIG"" | qsub -P jdhotopp-lab -l mem_free=5G -N novoplasty -wd "$OUTPUT_DIR"
```