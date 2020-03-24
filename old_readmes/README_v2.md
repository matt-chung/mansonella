# mansonella

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
  - [Software](#software)
  - [Directories](#directories)
  - [Create directories](#create-directories)
- [Assemble mitochondrial genome of Mansonella using 4 different samples](#assemble-mitochondrial-genome-of-mansonella-using-4-different-samples)
  - [Download Mansonella ozzardi mitochondra reference](#download-mansonella-ozzardi-mitochondra-reference)
  - [Assemble Mansonella mitochondrial genome with NOVOPlasty](#assemble-mansonella-mitochondrial-genome-with-novoplasty)
    - [Create contig files for NOVOPlasty assembly](#create-contig-files-for-novoplasty-assembly)
    - [Assemble Mansonella mitochondrial genome using short read FASTQs with NOVOPlasty](#assemble-mansonella-mitochondrial-genome-using-short-read-fastqs-with-novoplasty)
    - [Check if NOVOPlasty can assemble a complete Mansonella mitochondria without a reference sequence input](#check-if-novoplasty-can-assemble-a-complete-mansonella-mitochondria-without-a-reference-sequence-input)
    - [Identify proper T6 mitochondria assembly](#identify-proper-t6-mitochondria-assembly)
  - [Compare sequence identity between the the M. ozzardi mitogenome and the T6, T8 assemblies using BLASTN](#compare-sequence-identity-between-the-the-m-ozzardi-mitogenome-and-the-t6-t8-assemblies-using-blastn)
- [Create improved assembly of M. perstans mitochondria by combining long read data](#create-improved-assembly-of-m-perstans-mitochondria-by-combining-long-read-data)
  - [Identify T8 long and short reads that map to T8 NOVOPlasty assembly](#identify-t8-long-and-short-reads-that-map-to-t8-novoplasty-assembly)
    - [Map T8 short reads to T8 assembly](#map-t8-short-reads-to-t8-assembly)
    - [Create a subset FASTQ containing T8 short reads that mapped to T8 assembly](#create-a-subset-fastq-containing-t8-short-reads-that-mapped-to-t8-assembly)
    - [Map T8 long reads to T8 assembly](#map-t8-long-reads-to-t8-assembly)
    - [Create a subset FASTQ containing T8 short reads that mapped to T8 assembly](#create-a-subset-fastq-containing-t8-short-reads-that-mapped-to-t8-assembly-1)
  - [Assemble T8 mitochondria using subset short and long read FASTQs that contain reads that mapped to the T8 NOVOPlasty assembly](#assemble-t8-mitochondria-using-subset-short-and-long-read-fastqs-that-contain-reads-that-mapped-to-the-t8-novoplasty-assembly)
  - [T8](#t8)
  - [T6](#t6)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
PERL_BIN_DIR=/usr/local/packages/perl-5.24.0/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

BWA_BIN_DIR=/usr/local/packages/bwa-0.7.17/bin
MINIMAP2_BIN_DIR=/usr/local/packages/minimap2-2.17/bin
MUMMER_BIN_DIR=/usr/local/packages/mummer-3.23
NOVOPLASTY_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/NOVOPlasty_v3.8.1
PILON_BIN_DIR=/usr/local/packages/pilon-1.22
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
UNICYCLER_BIN_DIR=/usr/local/packages/python-3.5.2/bin
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
mkdir -p "$WORKING_DIR"/assemblies/unicycler/assembly1/P359_3

mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T6
mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T7_2
mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/T8
mkdir -p "$WORKING_DIR"/assemblies/novoplasty/assembly1/P359_3

mkdir -p "$WORKING_DIR"/nucmer/
mkdir -p "$WORKING_DIR"/final_assemblies/
```

# Assemble mitochondrial genome of Mansonella using 4 different samples

T6 = 5 larvae  
T8 = 3 larvae/adults  
T7.2 = 1 adult  
P359_3 = ???

## Download Mansonella ozzardi mitochondra reference
```{bash, eval = F}
wget 
```


## Assemble Mansonella mitochondrial genome with NOVOPlasty

### Create contig files for NOVOPlasty assembly

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
d
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

P359_3: "$WORKING_DIR"/assemblies/novoplasty/assembly1/P359_3/config.txt
```{bash, eval = F}
Project:
-----------------------
Project name          = P359_3
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
Forward reads         = /local/projects-t4/EMANS/P359_3_Manonella_Lib/ILLUMINA_DATA/EMANS_20180308_K00134_IL100099612_S5_L004_R1_trimmed.fastq.gz
Reverse reads         = /local/projects-t4/EMANS/P359_3_Manonella_Lib/ILLUMINA_DATA/EMANS_20180308_K00134_IL100099612_S5_L004_R2_trimmed.fastq.gz

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

### Assemble Mansonella mitochondrial genome using short read FASTQs with NOVOPlasty

Not enough coverage for mitochondrial assembly from T7_2 and P359_3

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

## P359_3
CONFIG="$WORKING_DIR"/assemblies/novoplasty/assembly1/P359_3/config.txt
OUTPUT_DIR="$WORKING_DIR"/assemblies/novoplasty/assembly1/P359_3/
``` 

##### Commands: 
```{bash, eval = F}
echo -e ""$PERL_BIN_DIR"/perl "$NOVOPLASTY_BIN_DIR"/NOVOPlasty3.8.1.pl -c "$CONFIG"" | qsub -P jdhotopp-lab -l mem_free=5G -N novoplasty -wd "$OUTPUT_DIR"
```

### Check if NOVOPlasty can assemble a complete Mansonella mitochondria without a reference sequence input

The assemblies differ by one ambigious base pair position.

##### Inputs:
```{bash, eval = F}
## T8
CONFIG="$WORKING_DIR"/assemblies/novoplasty/assembly1/T8_noref/config_noref.txt
OUTPUT_DIR="$WORKING_DIR"/assemblies/novoplasty/assembly1/T8_noref/
```

##### Commands: 
```{bash, eval = F}
echo -e ""$PERL_BIN_DIR"/perl "$NOVOPLASTY_BIN_DIR"/NOVOPlasty3.8.1.pl -c "$CONFIG"" | qsub -P jdhotopp-lab -l mem_free=5G -N novoplasty -wd "$OUTPUT_DIR"
```

No reference assembly stats:
```{bash, eval = F}
/home/jdhotopp/bin/residues.pl /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T8_noref/Circularized_assembly_1_T8.fasta
```
```{bash, eval = F}
Contig1 13621   13621   perG+C:25.8     T:7193  C:995   G:2520  A:2857  other:56
```

With reference assembly stats:
```{bash, eval = F}
/home/jdhotopp/bin/residues.pl /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T8/Circularized_assembly_1_T8.fasta
```
```{bash, eval = F}
Contig1 13622   13622   perG+C:25.8     T:7193  C:995   G:2520  A:2857  other:57
```

Difference in assemblies:

```{bash, eval = F}
tail -n1  /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T8/Circularized_assembly_1_T8.fasta
```
```{bash, eval = F}
AATTTTTTTATTTTATTAATGTTTTTTTGTTTTTTATATTTTTTGTTGTTTTTTTTTTTTT*GTTATGGTTAATTTTGTAATCATTTTATAATTTTTACTTTAGTAAAATTTGTTTTTGTGG
```
```{bash, eval = F}
tail -n1  /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T8_noref/Circularized_assembly_1_T8.fasta
```
```{bash, eval = F}
AATTTTTTTATTTTATTAATGTTTTTTTGTTTTTTATATTTTTTGTTGTTTTTTTTTTTTTGTTATGGTTAATTTTGTAATCATTTTATAATTTTTACTTTAGTAAAATTTGTTTTTGTGG
```

### Identify proper T6 mitochondria assembly

Only one circularized assembly was recovered from T8. There are four options (involving different combinations of four contigs) for T6:  

>Contig 01+03+04  
>Contig 01+02+05  
>Contig 01+02+04  
>Contig 01+03+05  

T8 look like they have a deletion at around 8.5 kb relative to the reference.

##### Inputs:
```{bash, eval = F}
REF_FNA=/local/projects-t3/EBMAL/mchung_dir/mansonella/references/KX822021.1.fna
T6_NOVOPLASTY_FNA=/local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T6/Contigs_1_T6.fasta
T8_NOVOPLASTY_FNA=/local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/novoplasty/assembly1/T8/Circularized_assembly_1_T8.fasta
```

##### Commands: 
```{bash, eval = F}
"$MUMMER_BIN_DIR"/nucmer "$REF_FNA" "$T6_NOVOPLASTY_FNA" -p "$WORKING_DIR"/nucmer/ref_v_novoT6
"$MUMMER_BIN_DIR"/nucmer "$REF_FNA" "$T8_NOVOPLASTY_FNA" -p "$WORKING_DIR"/nucmer/ref_v_novoT8
"$MUMMER_BIN_DIR"/nucmer "$T8_NOVOPLASTY_FNA" "$T6_NOVOPLASTY_FNA" -p "$WORKING_DIR"/nucmer/novoT8_v_novoT6
```

Reference v. NOVOPlasty T6 assembly:
```{bash, eval = F}
"$MUMMER_BIN_DIR"/mummerplot --layout --png "$WORKING_DIR"/nucmer/ref_v_novoT6.delta --prefix "$WORKING_DIR"/nucmer/ref_v_novoT6
```
![image](/images/ref_v_novoT6.png)

NOVOPlasty T8 assembly v. NOVOPlasty T6 assembly:
```{bash, eval = F}
"$MUMMER_BIN_DIR"/mummerplot --layout --png "$WORKING_DIR"/nucmer/novoT8_v_novoT6.delta --prefix "$WORKING_DIR"/nucmer/novoT8_v_novoT6
```
![image](/images/novoT8_v_novoT6.png)

Contig 01+02+05 (Option #2) looks like the correct T6 assembly.

## Compare sequence identity between the the M. ozzardi mitogenome and the T6, T8 assemblies using BLASTN

M. ozzardi genome used as query, T6 and T8 assemblies used as reference:

![image](/images/blast1.png)

T6 used as query, T8 used as reference:

![image](/images/blast2.png)

# Create improved assembly of M. perstans mitochondria by combining long read data

## Identify T8 long and short reads that map to T8 NOVOPlasty assembly

### Map T8 short reads to T8 assembly

##### Inputs:
```{bash, eval = F}
SEED_LENGTH=23
THREADS=16

## T8
OUTPUT_PREFIX=T8_v_novoT8
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
```

##### Commands:
```{bash, eval = F}
"$BWA_BIN_DIR"/bwa index "$REF_FNA"
echo -e ""$BWA_BIN_DIR"/bwa mem -t "$THREADS" -k "$SEED_LENGTH" "$REF_FNA" "$FASTQ1" "$FASTQ2" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N bwa -wd "$WORKING_DIR"/assemblies/
```
### Create a subset FASTQ containing T8 short reads that mapped to T8 assembly

##### Inputs:
```{bash, eval = F}
## T8
OUTPUT_PREFIX=T8_v_KX822021.1
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
BAM="$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F 4 "$BAM" | awk '{print $1}' | sort -n | uniq > "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".mappedreads.list
"$SEQTK_BIN_DIR"/seqtk subseq "$FASTQ1" "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".mappedreads.list | gzip > "$WORKING_DIR"/assemblies/"$(basename "$FASTQ1" | sed "s/[.]fastq.gz/.subset.fastq.gz/g")"
"$SEQTK_BIN_DIR"/seqtk subseq "$FASTQ2" "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".mappedreads.list | gzip > "$WORKING_DIR"/assemblies/"$(basename "$FASTQ2" | sed "s/[.]fastq.gz/.subset.fastq.gz/g")"
```

### Map T8 long reads to T8 assembly

##### Inputs:
```{bash, eval = F}
## T8
OUTPUT_PREFIX=T8_v_novoT8
REF_FNA="$WORKING_DIR"/final_assemblies/novoplasty_T8.fasta
FASTQ=/local/projects/RDMIN/SEQUENCE/20180806/fastq/emans_T8_mda_branched.fastq
```

##### Commands:
```{bash, eval = F}
"$MINIMAP2_BIN_DIR"/minimap2 -ax map-ont "$REF_FNA" "$FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".ont.bam -
```

### Create a subset FASTQ containing T8 short reads that mapped to T8 assembly

##### Inputs:
```{bash, eval = F}
## T8
OUTPUT_PREFIX=T8_v_novoT8
FASTQ=/local/projects/RDMIN/SEQUENCE/20180806/fastq/emans_T8_mda_branched.fastq
BAM="$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".ont.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools view -F 4 "$BAM" | awk '{print $1}' | sort -n | uniq > "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".ontmappedreads.list
"$SEQTK_BIN_DIR"/seqtk subseq "$FASTQ" "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".ontmappedreads.list | gzip > "$WORKING_DIR"/assemblies/"$(basename "$FASTQ" | sed "s/[.]fastq.gz/.subset.fastq.gz/g")"
```

## Assemble T8 mitochondria using subset short and long read FASTQs that contain reads that mapped to the T8 NOVOPlasty assembly

##### Inputs:
```{bash, eval = F}
THREADS=16
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
REF_GENE_FNA="$WORKING_DIR"/references/KX822021.1.gene.fna

## T8
OUTPUT_PREFIX=T8_v_novoT8
LONG_FASTQ="$WORKING_DIR"/assemblies/emans_T8_mda_branched.subset.fastq.gz
SHORT_FASTQ1="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.subset.fastq.gz
SHORT_FASTQ2="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.subset.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/assemblies/unicycler/assembly1/T8
```

##### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$UNICYCLER_BIN_DIR"/unicycler --mode bold --long "$LONG_FASTQ" --short1 "$SHORT_FASTQ1" --short2 "$SHORT_FASTQ2" -o "$OUTPUT_DIR" --pilon_path "$PILON_BIN_DIR"/pilon-1.22.jar -t "$THREADS"" | qsub -P jdhotopp-lab -N unicycler -wd "$OUTPUT_DIR" -q threaded.q -pe thread "$THREADS" -l mem_free=5G
```

qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N getorganelle -wd /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/getorganelle -b y /local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle/get_organelle_from_reads.py -1 "$SHORT_FASTQ1" -2 "$SHORT_FASTQ2" -s "$REF_FNA" --genes /local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle/GetOrganelleLib/LabelDatabase/animal_mt.fasta -R 10 -k 21,45,65,85,105 -F animal_mt -o /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/getorganelle -t "$THREADS"


qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N getorganelle -wd /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/getorganelle -b y /local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle/get_organelle_from_reads.py -1 "$SHORT_FASTQ1" -2 "$SHORT_FASTQ2" -s /local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle/GetOrganelleLib/SeedDatabase/animal_mt.fasta --genes /local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle/GetOrganelleLib/LabelDatabase/animal_mt.fasta -R 10 -k 21,45,65,85,105 -F animal_mt -o /local/projects-t3/EBMAL/mchung_dir/mansonella/assemblies/getorganelle -t "$THREADS"


## T8
OUTPUT_PREFIX=T8_v_getorganelleT8
REF_FNA="$WORKING_DIR"/assemblies/getorganelle/animal_mt.K105.complete.graph1.1.path_sequence.fasta
LONG_FASTQ="$WORKING_DIR"/assemblies/emans_T8_mda_branched.subset.fastq.gz
SHORT_FASTQ1="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.subset.fastq.gz
SHORT_FASTQ2="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.subset.fastq.gz

## T6
OUTPUT_PREFIX=T6_v_getorganelleT8
REF_FNA="$WORKING_DIR"/assemblies/getorganelle/animal_mt.K105.complete.graph1.1.path_sequence.fasta
SHORT_FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
SHORT_FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz


"$MINIMAP2_BIN_DIR"/minimap2 -ax map-ont "$REF_FNA" "$LONG_FASTQ" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".ont.bam -


/usr/local/packages/bwa-0.7.17/bin/bwa index "$REF_FNA"
 echo -e "/usr/local/packages/bwa-0.7.17/bin/bwa mem -t 4 -k 23 "$REF_FNA" "$SHORT_FASTQ1" "$SHORT_FASTQ2" | /usr/local/packages/samtools-1.3.1/bin/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".pe.bam -" | qsub -q threaded.q  -pe thread 4 -P jdhotopp-lab -l mem_free=5G -N bwa -wd "$(dirname $output_bam)"

