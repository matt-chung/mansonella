# Complete mitochondrial genome sequence of Mansonella perstans

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
  - [Software](#software)
  - [Directories](#directories)
- [Assemble M. ozzardi mitogenome](#assemble-m-ozzardi-mitogenome)
  - [Identify T6 and T8 long and short reads that map to M. ozzardi mitogenome](#identify-t6-and-t8-long-and-short-reads-that-map-to-m-ozzardi-mitogenome)
    - [Map T6 and T8 short reads to M. ozzardi mitogenome](#map-t6-and-t8-short-reads-to-m-ozzardi-mitogenome)
    - [Create a subset FASTQ containing T6 and T8 short reads that mapped to M. ozzardi mitogenome](#create-a-subset-fastq-containing-t6-and-t8-short-reads-that-mapped-to-m-ozzardi-mitogenome)
  - [Assemble T8 mitochondria using subset short read FASTQs that contain reads that mapped to M. ozzardi mitogenome](#assemble-t8-mitochondria-using-subset-short-read-fastqs-that-contain-reads-that-mapped-to-m-ozzardi-mitogenome)
  - [Rotate T8 assembly to M. ozzardi mitogenome](#rotate-t8-assembly-to-m-ozzardi-mitogenome)
    - [Run NUCMER on T8 assembly and M. ozzardi mitogenome](#run-nucmer-on-t8-assembly-and-m-ozzardi-mitogenome)
    - [Rotate T8 assembly from NUCMER coordinates](#rotate-t8-assembly-from-nucmer-coordinates)
- [Annotate the M. ozzardi and M. perstans mitogenome](#annotate-the-m-ozzardi-and-m-perstans-mitogenome)
  - [Use GeSeq to annotate the mitogenomes](#use-geseq-to-annotate-the-mitogenomes)
  - [Move GESeq output to file system](#move-geseq-output-to-file-system)
- [Compare synteny between the T8 and M. ozzardi mitogenomes using ACT](#compare-synteny-between-the-t8-and-m-ozzardi-mitogenomes-using-act)
  - [Create coords file readable in ACT](#create-coords-file-readable-in-act)
  - [Format GESeq GFF3 files for ACT](#format-geseq-gff3-files-for-act)
  - [Set file inputs for ACT](#set-file-inputs-for-act)
  - [Export ACT instance](#export-act-instance)
- [Construct a phylogenetic tree for the M. perstans and other nematode mitogenomes](#construct-a-phylogenetic-tree-for-the-m-perstans-and-other-nematode-mitogenomes)
  - [Prepare mitogenome reference files](#prepare-mitogenome-reference-files)
    - [Download nematode mitogenome references](#download-nematode-mitogenome-references)
    - [Copy T8 mitogenome assembly to references directory](#copy-t8-mitogenome-assembly-to-references-directory)
    - [Combine mitogenome FASTAs to a single multi-FASTA](#combine-mitogenome-fastas-to-a-single-multi-fasta)
  - [Align mitogenomes using MAFFT](#align-mitogenomes-using-mafft)
  - [Construct phylogenetic tree using IQ-TREE](#construct-phylogenetic-tree-using-iq-tree)
  - [Visualize phylogenetic tree using iTOL](#visualize-phylogenetic-tree-using-itol)
- [Assess core sequence identity between M. perstans mitogenome and 15 other filarial mitogenomes](#assess-core-sequence-identity-between-m-perstans-mitogenome-and-15-other-filarial-mitogenomes)
  - [Use mothur to filter MAFFT alignment to only contain positions present in all 16 mitogenomes](#use-mothur-to-filter-mafft-alignment-to-only-contain-positions-present-in-all-16-mitogenomes)
  - [Visualize core mitogenome sequence identity](#visualize-core-mitogenome-sequence-identity)
    - [Set R inputs](#set-r-inputs)
    - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
    - [Construct core sequence identity matrix](#construct-core-sequence-identity-matrix)
    - [Set heatmap rownames order and removes upper half of the triangle plot](#set-heatmap-rownames-order-and-removes-upper-half-of-the-triangle-plot)
    - [Plot sequence identity matrix as heatmap](#plot-sequence-identity-matrix-as-heatmap)
- [Identify SNPs and indels in the M. perstans mitogenome](#identify-snps-and-indels-in-the-m-perstans-mitogenome)
  - [Map T6 and T8 short reads to M. ozzardi mitogenome](#map-t6-and-t8-short-reads-to-m-ozzardi-mitogenome-1)
  - [Sort and index BAM files](#sort-and-index-bam-files)
  - [Generate MPILEUPs for BAM files](#generate-mpileups-for-bam-files)
  - [Count SNPs and indels in the MPILEUPs](#count-snps-and-indels-in-the-mpileups)
  - [Plot SNP and indel locations from T6 and T8 mapped against the T8 assembly](#plot-snp-and-indel-locations-from-t6-and-t8-mapped-against-the-t8-assembly)
    - [Set R inputs](#set-r-inputs-1)
    - [Load R functions and view sessionInfo](#load-r-functions-and-view-sessioninfo)
    - [Parse BASE files for plotting](#parse-base-files-for-plotting)
    - [Create depth plot with SNP and indel marks for each BASE file](#create-depth-plot-with-snp-and-indel-marks-for-each-base-file)
- [Construct final figures](#construct-final-figures)
  - [Figure 1](#figure-1)

<!-- /MarkdownTOC -->


# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software    
```{bash, eval = F}
PERL_BIN_DIR=/usr/local/packages/perl-5.24.0/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

ACT_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/artemis
BWA_BIN_DIR=/usr/local/packages/bwa-0.7.17/bin
EDIRECT_BIN_DIR=/usr/local/packages/edirect-9.19/
GETORGANELLE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle_v1.6.2e
IQTREE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/iqtree-1.6.2-Linux/bin
MAFFT_BIN_DIR=/usr/local/packages/mafft-7.427/bin
MOTHUR_BIN_DIR=/usr/local/packages/mothur-1.40.4
MUMMER_BIN_DIR=/usr/local/packages/mummer-3.23
NOVOPLASTY_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/NOVOPlasty_v3.8.1
PILEUP2BASE_BIN_DIR=/home/mattchung/scripts/pileup2base-master/
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin
```

## Directories

```{bash, eval = F}
REFERENCES_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella/references
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella
OUTPUT_DIR=/local/projects-t3/EBMAL/mchung_dir/mansonella/output
```

# Assemble M. ozzardi mitogenome

## Identify T6 and T8 long and short reads that map to M. ozzardi mitogenome

### Map T6 and T8 short reads to M. ozzardi mitogenome

##### Inputs:
```{bash, eval = F}
SEED_LENGTH=23
THREADS=16

## T6
OUTPUT_PREFIX=T6_v_KX822021.1
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz

## T8
OUTPUT_PREFIX=T8_v_KX822021.1
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
```

##### Commands:
```{bash, eval = F}
"$BWA_BIN_DIR"/bwa index "$REF_FNA"
echo -e ""$BWA_BIN_DIR"/bwa mem -t "$THREADS" -k "$SEED_LENGTH" "$REF_FNA" "$FASTQ1" "$FASTQ2" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N bwa -wd "$WORKING_DIR"/assemblies/
```

### Create a subset FASTQ containing T6 and T8 short reads that mapped to M. ozzardi mitogenome

##### Inputs:
```{bash, eval = F}
## T6
OUTPUT_PREFIX=T6_v_KX822021.1
FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz
BAM="$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".bam

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

## Assemble T8 mitochondria using subset short read FASTQs that contain reads that mapped to M. ozzardi mitogenome

##### Inputs:
```{bash, eval = F}
THREADS=4
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
REF_GENE_FNA="$WORKING_DIR"/references/KX822021.1.gene.fna

## T8
FASTQ1="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.subset.fastq.gz
FASTQ2="$WORKING_DIR"/assemblies/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.subset.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/assemblies/getorganelle/T8
```

##### Commands:
```{bash, eval = F}
"$GETORGANELLE_BIN_DIR"/get_organelle_from_reads.py -1 "$FASTQ1" -2 "$FASTQ2" -s "$REF_FNA" --genes "$REF_GENE_FNA" -R 10 -k 21,45,65,85,105 -F animal_mt -o "$OUTPUT_DIR" -t "$THREADS"
```

## Rotate T8 assembly to M. ozzardi mitogenome

##### Inputs
```{bash, eval = F}
OUTPUT_PREFIX=ref_v_getorganelleT8
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
GETORGANELLE_FNA="$WORKING_DIR"/assemblies/getorganelle/T8/animal_mt.K45.complete.graph1.1.path_sequence.fasta
```

### Run NUCMER on T8 assembly and M. ozzardi mitogenome

##### Commands
```{bash, eval = F}
"$MUMMER_BIN_DIR"/nucmer "$REF_FNA" "$GETORGANELLE_FNA" -p "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX"
"$MUMMER_BIN_DIR"/show-coords -bHT "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".delta
```

```{bash, eval = F}
1       2322    11327   13646   2322    2320    KX822021.1      164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+(circular)
2901    13611   613     11298   10711   10686   KX822021.1      164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+(circular)
```

### Rotate T8 assembly from NUCMER coordinates

T8 assembly is 13,647 bp long.

##### Commands
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools faidx "$GETORGANELLE_FNA" 164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+(circular):11327-13647 -o "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt1.fasta
"$SAMTOOLS_BIN_DIR"/samtools faidx "$GETORGANELLE_FNA" 164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+(circular):1-11298 -o "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt2.fasta

cat "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt1.fasta "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt2.fasta | grep -v ":1-11298" > "$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
rm "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt1.fasta
rm "$WORKING_DIR"/final_assemblies/"$OUTPUT_PREFIX"_pt2.fasta
```

# Annotate the M. ozzardi and M. perstans mitogenome

## Use GeSeq to annotate the mitogenomes

![image](/images/geseq_settings.png)

## Move GESeq output to file system

Unzipped job-results-20203684142 to "$WORKING_DIR"/geseq/job-results-20203684142

# Compare synteny between the T8 and M. ozzardi mitogenomes using ACT

## Create coords file readable in ACT

##### Inputs
```{bash, eval = F}
OUTPUT_PREFIX=ref_v_getorganelleT8
REF_FNA="$WORKING_DIR"/references/KX822021.1.fna
QUERY_FNA="$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
```

##### Commands
```{bash, eval = F}
"$MUMMER_BIN_DIR"/nucmer "$REF_FNA" "$QUERY_FNA" -p "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX"
"$MUMMER_BIN_DIR"/show-coords "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".delta > "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".coords
perl ~/scripts/nucmer_coords2ACT_galaxy.pl "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".coords "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".tab
```

## Format GESeq GFF3 files for ACT

##### Inputs
```{bash, eval = F}
## T8
GFF3="$WORKING_DIR"/geseq/job-results-20203684142/GeSeqJob-20200305-23934_164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+%28circular%29:11327-13647_GFF3.gff3

## M. ozzardi
GFF3="$WORKING_DIR"/geseq/job-results-20203684142/GeSeqJob-20200305-23934_KX822021.1-UNVERIFIED:-Mansonella-ozzardi-mitochondrion-sequence_GFF3.gff3
```

##### Commands
```{bash, eval = F}
awk '$2 == "ARWEN_v1.2.3" || $3 != "tRNA" {print $0}' "$GFF3" | awk '$3 == "CDS" || $3 == "rRNA" || $3 == "tRNA" {print $0}' > "$(echo "$GFF3" | sed "s/[.]gff3$/.act.gff3/g")"
```

## Set file inputs for ACT

##### Commands
```{bash, eval = F}
"$ARTEMIS_BIN_DIR"/act
```

```{bash, eval = F}
Sequence file 1: "$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
Comparison file 1: "$WORKING_DIR"/nucmer/"$OUTPUT_PREFIX".tab
Sequence file 2: "$WORKING_DIR"/references/KX822021.1.fna

Annotation file 1: "$WORKING_DIR"/geseq/job-results-20203684142/GeSeqJob-20200305-23934_164174_164176_163722_163582_163580_164124_156882_163258_164170_163684_164128_163726_163640_161910_162854_164116_158062_164154_156920_163106_164074_163638+%28circular%29:11327-13647_GFF3.act.gff3
Annotation file 2: "$WORKING_DIR"/geseq/job-results-20203684142/GeSeqJob-20200305-23934_KX822021.1-UNVERIFIED:-Mansonella-ozzardi-mitochondrion-sequence_GFF3.act.gff3
```

## Export ACT instance

![image](/images/act.png)

# Construct a phylogenetic tree for the M. perstans and other nematode mitogenomes

## Prepare mitogenome reference files

### Download nematode mitogenome references

##### Commands
```{bash, eval = F}
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id HM773029 > "$REFERENCES_DIR"/Chandlerella_quiscali_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id KX822021 > "$REFERENCES_DIR"/Mansonella_ozzardi_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id JN367461 > "$REFERENCES_DIR"/Wuchereria_bancrofti_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id AP017680 > "$REFERENCES_DIR"/Brugia_pahangi_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id AP017686 > "$REFERENCES_DIR"/Brugia_timori_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id AF538716 > "$REFERENCES_DIR"/Brugia_malayi_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id HQ186250 > "$REFERENCES_DIR"/Loa_loa_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id KX265048 > "$REFERENCES_DIR"/Dirofilaria_repens_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id AJ537512 > "$REFERENCES_DIR"/Dirofilaria_immitis_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id KT599912 > "$REFERENCES_DIR"/Onchocerca_volvulus_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id KX181289 > "$REFERENCES_DIR"/Onchocerca_ochengi_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id HQ214004 > "$REFERENCES_DIR"/Onchocerca_flexuosa_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id HQ186249 > "$REFERENCES_DIR"/Acanthocheilonema_viteae_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id GU138699 > "$REFERENCES_DIR"/Setaria_digitata_mitogenome.fna
"$EDIRECT_BIN_DIR"/efetch -format fasta -db nuccore -id MH937750 > "$REFERENCES_DIR"/Setaria_labiatopapillosa_mitogenome.fna
```

### Copy T8 mitogenome assembly to references directory

##### Commands
```{bash, eval = F}
cp "$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta "$REFERENCES_DIR"/Mansonella_perstans_mitogenome.fna
```

### Combine mitogenome FASTAs to a single multi-FASTA

##### Commands
```{bash, eval = F}
cat "$REFERENCES_DIR"/*mitogenome.fna > "$REFERENCES_DIR"/mitogenome_combined.fna
```

## Align mitogenomes using MAFFT

##### Inputs
```{bash, eval = F}
FNA="$REFERENCES_DIR"/mitogenome_combined.fna
```

##### Commands
```{bash, eval = F}
"$MAFFT_BIN_DIR"/mafft "$FNA" > "$(echo "$FNA" | sed "s/[.]fna*/.aligned.fna/g")"
```

## Construct phylogenetic tree using IQ-TREE
##### Inputs
```{bash, eval = F}
THREADS=4
MSA_FNA="$REFERENCES_DIR"/mitogenome_combined.aligned.fna
```
##### Commands
```{bash, eval = F}
"$IQTREE_BIN_DIR"/iqtree -s "$MSA_FNA" -nt "$THREADS" -bb 1000 -redo
```

## Visualize phylogenetic tree using iTOL

Site: itol.embl.de  

TREE_FILE: "$REFERENCES_DIR"/mitogenome_combined.aligned.fna.treefile

![image](/images/itol.png)

# Assess core sequence identity between M. perstans mitogenome and 15 other filarial mitogenomes

## Use mothur to filter MAFFT alignment to only contain positions present in all 16 mitogenomes

##### Inputs
```{bash, eval = F}
MSA_FNA="$REFERENCES_DIR"/mitogenome_combined.aligned.fna
```

##### Commands
```{bash, eval = F}
"$MOTHUR_BIN_DIR"/mothur "#filter.seqs(fasta="$MSA_FNA", vertical=F, trump=-)"
"$MOTHUR_BIN_DIR"/mothur "#filter.seqs(fasta="$(echo "$MSA_FNA" | sed "s/[.]fna$/filter.fasta/g")", vertical=F, trump=.)"
mv "$(echo "$MSA_FNA" | sed "s/[.]fna$/filter.fasta/g")" "$(echo "$MSA_FNA" | sed "s/[.]fna$/core.fna/g")"
```

## Visualize core mitogenome sequence identity

### Set R inputs
```{R}
CORE_FNA.PATH <- "Z:EBMAL/mchung_dir/mansonella/references/mitogenome_combined.aligned.core.fna"
OUTPUT.DIR <- "Z:EBMAL/mchung_dir/mansonella/plots"
```

### Load R packages and view sessionInfo
```{r}
library(Biostrings)
library(ggplot2)
library(pvclust)
library(reshape)

sessionInfo()
```

```{r, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape_0.8.8       pvclust_2.0-0       ggplot2_3.2.0       Biostrings_2.50.2   XVector_0.22.0      IRanges_2.16.0     
[7] S4Vectors_0.20.1    BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       rstudioapi_0.10  knitr_1.23       magrittr_1.5     zlibbioc_1.28.0  tidyselect_0.2.5 munsell_0.5.0   
 [8] colorspace_1.4-1 R6_2.4.0         rlang_0.4.0      plyr_1.8.4       dplyr_0.8.3      tools_3.5.0      grid_3.5.0      
[15] gtable_0.3.0     xfun_0.8         withr_2.1.2      yaml_2.2.0       lazyeval_0.2.2   assertthat_0.2.1 tibble_2.1.3    
[22] crayon_1.3.4     purrr_0.3.2      glue_1.3.1       compiler_3.5.0   pillar_1.4.2     scales_1.0.0     pkgconfig_2.0.2 
```

### Construct core sequence identity matrix
```{r}
concat.alignment <- readDNAMultipleAlignment(CORE_FNA.PATH)

sim <- as.data.frame(matrix(nrow = length(concat.alignment@unmasked),
                            ncol = length(concat.alignment@unmasked)))
rownames(sim) <- names(concat.alignment@unmasked)
rownames(sim) <- gsub("_.*","",gsub(" ","_",rownames(sim)))

colnames(sim) <- rownames(sim)

for(i in 1:nrow(sim)){
  j <- 1
  while(j <= i){
    palign <- PairwiseAlignments(c(concat.alignment@unmasked[i],concat.alignment@unmasked[j]))
    sim[i,j] <- pid(palign,type="PID4")
    j <- j+1
  }
}

for(i in 1:nrow(sim)){
  sim[i,] <- t(sim[,i])
}
```

### Set heatmap rownames order and removes upper half of the triangle plot
```{r}
levels <- pvclust(sim, nboot=100)
levels <- rownames(sim)[levels$hclust$order]

sim <- sim[match(levels,rownames(sim)),match(levels,colnames(sim))]

sim[lower.tri(sim)] <- NA
```

### Plot sequence identity matrix as heatmap
```{r, fig.width=7,fig.height=5}
plot.df <- as.data.frame(cbind(rownames(sim),
                               sim))
names(plot.df)[1] <- "temp"
plot.df <- melt(plot.df,id.vars="temp", na.rm=T)
colnames(plot.df) <- c("Var1", "Var2", "value")
plot.df$Var1 <- factor(plot.df$Var1, levels=levels)
plot.df$Var2 <- factor(plot.df$Var2, levels=rev(levels))

cgasi.hm <- ggplot(data = plot.df, aes(Var1, Var2, fill = value))+
  geom_tile(color = "black")+
  geom_text(aes(label = ifelse(value>=85, round(value, 1), "")), size = 3)+
  scale_fill_gradientn(colors = colorRampPalette(c("blue","darkturquoise","darkgreen","green","yellow","red"))(11),
                       limits=c(0,100))+
  theme_minimal()+
  labs(x="",y="",fill="Core Sequence\nIdentity")+
  #guides(fill = F)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1,size=8),
        axis.text.y = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(paste0(OUTPUT.DIR,"/Fig1b_raw.pdf"),
    height=5,
    width=7)
print(cgasi.hm)
dev.off()

png(paste0(OUTPUT.DIR,"/Fig1b_raw.png"),
    height=5,
    width=7,
    units = "in",res=300)
print(cgasi.hm)
dev.off()

print(cgasi.hm)
```

# Identify SNPs and indels in the M. perstans mitogenome

## Map T6 and T8 short reads to M. ozzardi mitogenome

##### Inputs:
```{bash, eval = F}
SEED_LENGTH=23
THREADS=16

## T6
OUTPUT_PREFIX=T6_v_getorganelleT8
REF_FNA="$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
FASTQ1=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T6/ILLUMINA_DATA/EMANS_20181012_K00134_IL100106386_S1_L001_R2_trimmed.fastq.gz

## T8
OUTPUT_PREFIX=T8_v_getorganelleT8
REF_FNA="$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
FASTQ1=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R1_trimmed.fastq.gz
FASTQ2=/local/projects/EMANS/T8/ILLUMINA_DATA/EMANS_20180815_K00134_IL100106041_S21_L005_R2_trimmed.fastq.gz
```

##### Commands:
```{bash, eval = F}
"$BWA_BIN_DIR"/bwa index "$REF_FNA"
echo -e ""$BWA_BIN_DIR"/bwa mem -t "$THREADS" -k "$SEED_LENGTH" "$REF_FNA" "$FASTQ1" "$FASTQ2" | "$SAMTOOLS_BIN_DIR"/samtools view -bho "$WORKING_DIR"/assemblies/"$OUTPUT_PREFIX".bam -" | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=5G -N bwa -wd "$WORKING_DIR"/assemblies/
```

## Sort and index BAM files

##### Input Sets:
```{bash, eval = F}
THREADS=4

## T6
BAM="$WORKING_DIR"/assemblies/T6_v_getorganelleT8.bam

## T8
BAM="$WORKING_DIR"/assemblies/T8_v_getorganelleT8.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools sort -@ "$THREADS" -o "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/g")" "$BAM"
"$SAMTOOLS_BIN_DIR"/samtools index -@ "$THREADS" "$(echo "$BAM" | sed "s/[.]bam$/.sortedbyposition.bam/g")"
```

## Generate MPILEUPs for BAM files

##### Input Sets:
```{bash, eval = F}
REF_FNA="$WORKING_DIR"/final_assemblies/getorganelle_T8.fasta
THREADS=4

## T6
BAM="$WORKING_DIR"/assemblies/T6_v_getorganelleT8.sortedbyposition.bam

## T8
BAM="$WORKING_DIR"/assemblies/T8_v_getorganelleT8.sortedbyposition.bam
```

##### Commands:
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools mpileup -f -aa -d 1000000 -o "$(echo "$BAM" | sed "s/[.]bam$/.mpileup/g")" "$BAM"
```

## Count SNPs and indels in the MPILEUPs
##### Input Sets:
```{bash, eval = F}
## T6
MPILEUP="$WORKING_DIR"/assemblies/T6_v_getorganelleT8.sortedbyposition.mpileup

## T8
MPILEUP="$WORKING_DIR"/assemblies/T8_v_getorganelleT8.sortedbyposition.mpileup
```

##### Commands:
```{bash, eval = F}
"$PILEUP2BASE_BIN_DIR"/pileup2baseindel_no_strand.pl "$MPILEUP" 0 "$(echo "$MPILEUP" | sed "s/[.]mpileup$/.base/g")"
```

## Plot SNP and indel locations from T6 and T8 mapped against the T8 assembly

### Set R inputs
```{R}
BASE.DIR <- "Z:/EBMAL/mchung_dir/mansonella/assemblies/"
OUTPUT.DIR <- "Z:EBMAL/mchung_dir/mansonella/plots"
```

### Load R functions and view sessionInfo
```{R}
library(cowplot)
library(ggplot2)
library(reshape)
library(scales)

sessionInfo()
```

```{R, eval = F}
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] scales_1.0.0  reshape_0.8.8 ggplot2_3.2.0 cowplot_1.0.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       rstudioapi_0.10  knitr_1.23       magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0    colorspace_1.4-1
 [8] R6_2.4.0         rlang_0.4.0      plyr_1.8.4       dplyr_0.8.3      tools_3.5.0      grid_3.5.0       gtable_0.3.0    
[15] xfun_0.8         withr_2.1.2      yaml_2.2.0       lazyeval_0.2.2   assertthat_0.2.1 tibble_2.1.3     crayon_1.3.4    
[22] purrr_0.3.2      glue_1.3.1       compiler_3.5.0   pillar_1.4.2     pkgconfig_2.0.2
```

### Parse BASE files for plotting
```{R}
mpileup_base.files <- list.files(BASE.DIR, pattern="*[.]base$", full.names = T)
varcutoff <- 0.25

plot.df <- list()

for(i in 1:length(mpileup_base.files)){
  mpileup_base_sum <- read.delim(mpileup_base.files[i],stringsAsFactors = F)
  mpileup_base_sum$Insertion <- unlist(apply(mpileup_base_sum,1,function(x){
    if(!is.na(x[8])){
      return(sum(as.numeric(as.character(gsub(":.*","",unlist(strsplit(as.character(x[8]),split="[|]")))))))
    }else{
      return(0)
    }
  }))
  mpileup_base_sum$Deletion <- unlist(apply(mpileup_base_sum,1,function(x){
    if(!is.na(x[9])){
      return(sum(as.numeric(as.character(gsub(":.*","",unlist(strsplit(as.character(x[9]),split="[|]")))))))
    }else{
      return(0)
    }
  }))
  
  mpileup_base_sum$depth <- rowSums(mpileup_base_sum[,4:9])
  
  mpileup_base_sum$snppos <- apply(mpileup_base_sum,1,function(x){
    snp_cutoff <- as.numeric(as.character(x[10]))*varcutoff
    ifelse(length(which(as.numeric(as.character(x[4:7])) > snp_cutoff)) > 1,return(T),return(F))
  })
  
  mpileup_base_sum$inspos <- apply(mpileup_base_sum,1,function(x){
    indel_cutoff <- as.numeric(as.character(x[10]))*varcutoff
    ifelse(length(which(as.numeric(as.character(x[8])) > indel_cutoff)) > 0,return(T),return(F))
  })
  mpileup_base_sum$delpos <- apply(mpileup_base_sum,1,function(x){
    indel_cutoff <- as.numeric(as.character(x[10]))*varcutoff
    ifelse(length(which(as.numeric(as.character(x[9])) > indel_cutoff)) > 0,return(T),return(F))
  })
  
  plot.df[[i]] <- mpileup_base_sum
}
```

### Create depth plot with SNP and indel marks for each BASE file
```{R,fig.width=20,fig.height=3}
plot.list <- list()
for(i in 1:length(plot.df)){
  snppos <- plot.df[[i]][which(plot.df[[i]]$snppos == T),c(2,4,5,6,7,10)]
  snppos <- melt(snppos,id.vars = c("loc"),measure.vars = c("A","T","G","C"))
  inspos <- plot.df[[i]][which(plot.df[[i]]$inspos == T),c(2,4,5,6,7,10)]
  delpos <- plot.df[[i]][which(plot.df[[i]]$delpos == T),c(2,4,5,6,7,10)]

  plot.list[[i]] <- ggplot()+
    geom_ribbon(mapping=aes_string(x=plot.df[[i]]$loc,ymin=0,ymax=plot.df[[i]]$depth),fill = "grey")+
    geom_bar(mapping=aes_string(x=snppos$loc,y=snppos$value,fill=snppos$variable),width=5,stat="identity")+
    scale_fill_manual(values=c("green","red","orange","blue"))+
    scale_x_continuous(expand=c(0,0),label=comma)+
    scale_y_continuous(expand=c(0,0),label=comma)+
    labs(x="position",y="sequencing depth")+
    guides(fill = F)+
    theme_bw()
  
  if(nrow(inspos) > 0){
    plot.list[[i]] <- plot.list[[i]]+
      geom_point(mapping=aes_string(x=inspos$loc,y=0),shape=15)
  }
  if(nrow(delpos) > 0){
    plot.list[[i]] <- plot.list[[i]]+
      geom_point(mapping=aes_string(x=delpos$loc,y=0),shape=15,color="red")
  }
}

pdf(paste0(OUTPUT.DIR,"/Fig1b_pt1.pdf"),
    height=3,
    width=20)
plot_grid(plotlist = plot.list,nrow=2,ncol=1)
dev.off()

png(paste0(OUTPUT.DIR,"/Fig1b_pt1.png"),
    height=3,
    width=20,
    units = "in",res=300)
plot_grid(plotlist = plot.list,nrow=2,ncol=1)
dev.off()

plot_grid(plotlist = plot.list,nrow=2,ncol=1)
```

![image](/images/Fig1b_pt1.png)

# Construct final figures

## Figure 1

![image](/images/Fig1.png)
