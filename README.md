# RNASeq and Call Variants

* Requirements

- STAR v2.7.5c
- RSEM v1.3.1
- human genome GRCh37.75


# Make References

```bash
# get reference genome fasta
wget -c ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

# get refernece genome gtf
wget -c ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

gunzip -v Homo_sapiens.GRCh37.75.gtf.gz
gunzip -v Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

GTF="Homo_sapiens.GRCh37.75.gtf"
FASTA="Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"

mkdir -p rsem

# rsem: prepare reference
/scratch/apps/RSEM/rsem-prepare-reference --gtf $GTF \
        --bowtie2 \
        $FASTA \
        rsem/rsem

# caminho do executavel do programa STAR
GENOMEDIR="/scratch/refs/release-75/"
STAR="/scratch/apps/STAR/bin/Linux_x86_64/STAR"
RSEM="/scratch/apps/RSEM"

# 1/4 - Generating genome indexes
$STAR --runThreadN 6 \
       --runMode genomeGenerate \
       --genomeDir $GENOMEDIR \
       --genomeFastaFiles $FASTA \
       --sjdbGTFfile $GTF \
       --sjdbOverhang 149
```

# STAR: mapping and counting

```bash

### STAR

# paths
GENOMEDIR="/scratch/refs/release-75"
STAR="/scratch/apps/STAR/bin/Linux_x86_64/STAR"
RSEM="/scratch/apps/RSEM"
FASTA="$GENOMEDIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa "
GTF="$GENOMEDIR/Homo_sapiens.GRCh37.75.gtf"

# create a file with unique sample names
ls  -1d fastq/*/*.fastq.gz |  grep  -v  "Undetermined\|PK01" | sed  -e  "s/_L00[0-4]_R[12]_001.fastq.gz//g" | sort  -Vu > FileSample

# create directory for STAR output: RNASEQ_data
mkdir -p  RNASEQ_data
mkdir -p ~/RNASEQ_data/


# run STAR mapping and counting
for f in  `cat FileSample`;
do
        ls $f/*/*R1*.fastq.gz | tr '\n' ',' | sed -s "s/\,$/ /g" > R1
        ls $f/*/*R2*.fastq.gz | tr '\n' ',' | sed -s "s/\,$/ /g" > R2

        fastqR1=$(cat R1)
        fastqR2=$(cat R2)

        mkdir  "RNASEQ_data/$(basename $f)";
        $STAR --genomeDir $GENOMEDIR \
        --readFilesCommand zcat \
        --readFilesIn $fastqR1 $fastqR2 \
        --limitBAMsortRAM 8000000000 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --twopassMode Basic \
        --outFilterMultimapNmax  1 \
        --quantMode TranscriptomeSAM \
        --runThreadN 1 \
        --outFileNamePrefix "RNASEQ_data/$(basename $f)/";

        rm -f R1 R2

        mkdir "RNASEQ_data/rsem.$(basename $f)";
        $RSEM/rsem-calculate-expression --bam --no-bam-output -p 5 \
        --paired-end \
        --forward-prob 0 \
        RNASEQ_data/$(basename $f)/Aligned.toTranscriptome.out.bam \
        $GENOMEDIR/rsem/rsem \
        RNASEQ_data/rsem.$(basename $f)/rsem;

        #rm -rf "RNASEQ_data/$(basename $f)";
        #rm -f R1 R2

        rsync --progress -r "RNASEQ_data/$(basename $f)" ~/RNASEQ_data/$(basename $f)

        rm -rf "RNASEQ_data/$(basename $f)";
done

exit

# by genes
cd RNASEQ_data

mkdir gene-level
cd gene-level
ls -d1 ../rsem.* | gawk '{print("ln -s",$1"/rsem.genes.results",gensub("../rsem.","","g",$1))}' | sh
cd ..
time R --file=../run.merge.files.R --args gene-level 5 gene-level-5

# by isoforms
mkdir transcript-level
cd transcript-level
ls -d1 ../rsem.* | gawk '{print("ln -s",$1"/rsem.genes.results",gensub("../rsem.","","g",$1))}' | sh
cd ..
time R --file=../run.merge.files.R --args transcript-level 5 transcript-level
```
