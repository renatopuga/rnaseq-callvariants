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



