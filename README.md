# RNASeq and Call Variants
> The data will be available after publication.

## Machine 

- Ubuntu Desktop 20.04 LTS
- 6 cpus
- 32GB RAM Memory
- 256GB SSD (/scratch) and 1TB HDD

## Requirements

- Docker version 19.03.8
- STAR v2.7.5c
- RSEM v1.3.1
- human genome GRCh37.75
- samtools 1.10
- vcf-merge
- VEP ensemble v102

## Sample Information

- TruSeq stranded mRNA 150 cycles 
- NextSeq 550

| Lineage | Phenotype |
| ------- | --------- |
| PK03    | Controle  |
| PK05    | Controle  |
| PK06    | Controle  |
| CQ16    | Controle  |
| PK01    | PSP       |
| PK02    | PSP       |
| PK04    | PSP       |
| PK07    | PSP       |
| PK08    | PSP       |
| PK09    | PSP       |

## References: Download and Compile

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

## Mapping and counting (STAR and RSEM)

```bash

### STAR

# paths
GENOMEDIR="/scratch/refs/release-75"
STAR="/scratch/apps/STAR/bin/Linux_x86_64/STAR"
RSEM="/scratch/apps/RSEM"
FASTA="$GENOMEDIR/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa "
GTF="$GENOMEDIR/Homo_sapiens.GRCh37.75.gtf"

# create a file with unique sample names
ls  -1d fastq/*/*.fastq.gz |  grep  -v  "Undetermined" | sed  -e  "s/_L00[0-4]_R[12]_001.fastq.gz//g" | sort  -Vu > FileSample

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

## EdgeR - Differential Gene Expression

```bash
# Fonte: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# carregar a biblioteca edgeR
library("edgeR")

# table with gene count
# First Column: Symbol
# The remaining columns are the gene count
TabelaCount <- "merge-table-STAR-gene-level-5-50x.csv"

# loading the table (the first column is called Symbol)
x <- read.delim(TabelaCount,row.names="symbol")

# groups: 1 and 2
# The first 4 samples group 1
# 5 afterwards is group 2
group <- factor(c(1,1,1,1,2,2,2,2,2))
y <- DGEList(counts=x,group=group)


# We filter out lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# TMM normalization and display the normalization factors
y <- calcNormFactors(y)

png("plotMDS.png")
plotMDS(y)
dev.off()

design <- model.matrix(~group)
y <- estimateDisp(y,design)

png("plotBCV.png")
plotBCV(y)
dev.off()

#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)

png("plotMD.png")
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()

write.table(lrt$table,"lrt.DGElist.csv",quote=FALSE, sep="\t")
```

## GATK4 - Call Variants
> by https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels

```bash
# scratch
genome="/scratch/refs/release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
dbsnp="/scratch/bucket/dbsnp_138.b37.vcf.gz"
intervals="/scratch/refs/release-75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.interval_list"

tmp="/data/tmp"

sample=$1
LB="RNASeq"
PL="illumina"
PU="NextSeq500"

# listar amostras
bamList=$(ls -1 */*/*.bam)
for bam in $bamList
do
   out=$(echo $bam | sed "s/\/.*//")
 
# run MarkDuplicates
docker run  --user "$(id -u):$(id -g)" -it --rm -v /tmp:/tmp -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk --java-options "-Djava.io.tmpdir=${tmp}  -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" MarkDuplicates \
	--TMP_DIR $tmp \
	-I /data/$bam -O /data/$out.sorted.dup.bam \
       -M /data/$out.sorted.dup_metrics \
       --VALIDATION_STRINGENCY SILENT \
	--CREATE_INDEX true \


# run SplitNCigarReads
docker run  --user "$(id -u):$(id -g)" -it --rm -v /tmp:/tmp -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk SplitNCigarReads \
	-R $genome \
	-L $intervals_genes \
	-I /data/$out.sorted.dup.bam \
	-O /data/$out.sorted.dup.split.bam

# run AddOrReplaceGroups
docker run  --user "$(id -u):$(id -g)" -it --rm -v /tmp:/tmp -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk AddOrReplaceReadGroups \
	-I /data/$out.sorted.dup.split.bam \
	-O /data/$out.sorted.dup.split.gp.bam \
	-ID $out \
	-LB $LB \
	-PL $PL\
	-PU $PU \
	-SM $out

# run index ...gp.bam
samtools index $out.sorted.dup.split.gp.bam

# run BaseRecalibrator
docker run --user "$(id -u):$(id -g)" -it --rm -v /scratch/:/scratch -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator \
	-L $intervals_genes \
	-R $genome \
	-I /data/$out.sorted.dup.split.gp.bam \
	--known-sites $dbsnp \
	-O /data/$out.sorted.dup.split.gp.recal.table

# run ApplyBQSR
docker run --user "$(id -u):$(id -g)" -it --rm -v /scratch:/scratch -v $(pwd):/data broadinstitute/gatk gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR \
	-R $genome \
	-I /data/$out.sorted.dup.split.gp.bam \
	-bqsr /data/$out.sorted.dup.split.gp.recal.table \
	-L $intervals_genes --create-output-bam-index true --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
	-O /data/$out.sorted.dup.split.gp.recal.bam

# run HaplotypeCaller
docker run  --user "$(id -u):$(id -g)" -it --rm -v /scratch:/scratch  -v $(pwd):/data broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" HaplotypeCaller \
        -R $genome \
	-I /data/$out.sorted.dup.split.gp.recal.bam \
        -L $intervals_genes \
	-dont-use-soft-clipped-bases \
	--standard-min-confidence-threshold-for-calling 20 \
        --native-pair-hmm-threads 8 \
        -O /data/$out.g.vcf.gz \
	--dbsnp $dbsnp

# run VariantFiltration
docker run --user "$(id -u):$(id -g)" -it --rm -v /scratch:/scratch -v $(pwd):/data/ broadinstitute/gatk:4.1.4.1 gatk VariantFiltration \
        --R $genome \
        --V /data/$out.g.vcf.gz \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O /data/$out.vf.vcf.gz

done

# run list vcfs
vcfs=$(ls -1 *.vf.vcf.gz)
for i in $vcfs
do
	vcfList=${vcfList}" ${i}"
done

# run vcf-merge
vcf-merge $vcfList > samples.vcf

# create vep_data output
mkdir vep_data
chmod a+w vep_data

# run vep annotation
docker run -it --rm  -v $(pwd):/data -v /scratch/:/scratch -v /vep/:/opt/vep/.vep  ensemblorg/ensembl-vep ./vep  \
	-i /data/samples.vcf  \
	-o /data/vep_data/samples.vep.tsv -v \
	--fork 4 --cache --distance 0 --use_given_ref \
	--pick --pick_allele --no_intergenic --exclude_predicted \
	--no_stats --offline --force_overwrite --everything --tab \
	--fasta $genome  \
	--individual all 
```
