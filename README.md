# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
PERL_BIN_DIR=/usr/local/packages/perl-5.24.0/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

BWA_BIN_DIR=/usr/local/packages/bwa-0.7.17/bin
GETORGANELLE_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/GetOrganelle_v1.6.2e
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

## Identify T8 long and short reads that map to M. ozzardi mitogenome

### Map T8 short reads to M. ozzardi mitogenome

##### Inputs:
```{bash, eval = F}
SEED_LENGTH=23
THREADS=16

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
### Create a subset FASTQ containing T8 short reads that mapped to M. ozzardi mitogenome

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




export LD_LIBRARY_PATH=/usr/local/packages/qt-5.8.0/lib/:"$LD_LIBRARY_PATH"


##### Commands:
```{bash, eval = F}
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$UNICYCLER_BIN_DIR"/unicycler --mode bold --long "$LONG_FASTQ" --short1 "$SHORT_FASTQ1" --short2 "$SHORT_FASTQ2" -o "$OUTPUT_DIR" --pilon_path "$PILON_BIN_DIR"/pilon-1.22.jar -t "$THREADS"" | qsub -P jdhotopp-lab -N unicycler -wd "$OUTPUT_DIR" -q threaded.q -pe thread "$THREADS" -l mem_free=5G
```