#!/bin/bash

bash code/header.sh

DRAFT_ASSEMBLY=phillippy_plus_cortex_method/sampleDataOxford/windows/polished.fa
MISEQ_READS="/Net/cycloid/data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R1_001_val_1.fq /Net/cycloid/data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R2_001_val_2.fq"
DIR=phillippy_plus_cortex_method/sampleDataOxford/5kb_windows/
PREFIX=sampleDataOxford
WINDOW_SIZE=5000
NUM_PROCS=20
k=31
label='stringent'

# 1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)
echo "1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)"
mkdir -p $DIR/windows_results
cp phillippy_plus_cortex_method/sampleDataOxford/polished.fa $DIR
#bwa index $DIR/polished.fa
#bwa mem $DIR/polished.fa $MISEQ_READS > $DIR/$PREFIX.sam
#samtools view -Sb $DIR/$PREFIX.sam > $DIR/$PREFIX.bam
#samtools sort -o $DIR/$PREFIX_sorted.bam $DIR/$PREFIX.bam
cp phillippy_plus_cortex_method/sampleDataOxford/windows/"$PREFIX"_sorted.bam $DIR/"$PREFIX"_sorted.bam
samtools index -b $DIR/"$PREFIX"_sorted.bam

# 2. Divide up into windows
echo "2. Divide up into $WINDOW_SIZE b windows"
#mkdir -p $DIR/windows_results
contigs=$(grep \> $DRAFT_ASSEMBLY |cut -c 2-)
for contig in $contigs; do
	contig_read=$(sed -n "/$contig/{n;p;}" $DRAFT_ASSEMBLY)
	len_contig=$(sed -n "/$contig/{n;p;}" $DRAFT_ASSEMBLY | wc -c)
	start_pos=0
	while [ $start_pos -lt $len_contig ]; do
		end_pos=$((start_pos + WINDOW_SIZE))
		mkdir -p $DIR/windows_results/$contig.$start_pos\-$end_pos/k$k$label
		contig_read_substring=${contig_read:$start_pos:$WINDOW_SIZE}
		echo -e ">$contig.$start_pos-$end_pos\n$contig_read_substring" > $DIR/windows_results/$contig.$start_pos-$end_pos/k$k$label/draft_assembly.$contig.$start_pos-$end_pos.fa
		samtools view  -b $DIR/$PREFIX"_sorted.bam" $contig:$start_pos-$end_pos -o $DIR/windows_results/$contig.$start_pos-$end_pos/k$k$label/$PREFIX.$contig.$start_pos-$end_pos.bam
		samtools bam2fq $DIR/windows_results/$contig.$start_pos-$end_pos/k$k$label/$PREFIX.$contig.$start_pos-$end_pos.bam > $DIR/windows_results/$contig.$start_pos-$end_pos/k$k$label/$PREFIX.$contig.$start_pos-$end_pos.fq
		let start_pos=$end_pos
	done
done

# 3. For each window, run_calls and correct
echo "3. For each window, run_calls and correct"
ls -d $PWD/$DIR/windows_results/*/k$k$label | parallel --gnu -j $NUM_PROCS 'bash code/cortex_correction.sh {} {}/draft_assembly*.fa '$WINDOW_SIZE' '$k' {}/*.fq &>>{}/log_cortex_correction_k'$k$label'.txt'

# 4. Combine corrected windows
echo "4. Combine corrected windows"
for contig in $contigs; do
        contig_read=$(sed -n "/$contig/{n;p;}" $DRAFT_ASSEMBLY)
        len_contig=$(sed -n "/$contig/{n;p;}" $DRAFT_ASSEMBLY | wc -c)
        start_pos=0
        while [ $start_pos -lt $len_contig ]; do
		if [ $start_pos -eq 0 ]; then
			echo ">$contig" &>> $DIR/corrected_assembly_k$k$label.fa
		fi
		end_pos=$((start_pos + WINDOW_SIZE))
		# want all but first line, and no newlines
		chunk_seq=$(tail -n +2 $DIR/windows_results/$contig.$start_pos-$end_pos/k$k$label/ref_bc_k$k.fa)
		#len_chunk=`expr length $chunk_seq`
		len_chunk=${#chunk_seq} 
		if [ "$len_chunk" -lt "$WINDOW_SIZE" ]; then
			echo "$start_pos seq is $len_chunk bps long, shorter than $WINDOW_SIZE"
		fi
		echo "$chunk_seq"  &>> $DIR/corrected_assembly_k$k$label.fa.tmp
		let start_pos=$end_pos
        done
done
cat $DIR/corrected_assembly_k$k$label.fa.tmp | tr -d '[:space:]' &>> $DIR/corrected_assembly_k$k$label.fa
rm $DIR/corrected_assembly_k$k$label.fa.tmp

bash code/footer.sh
