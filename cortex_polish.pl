#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $header = "bash code/header2.sh";
my $rheader = qx{$header};

#my $draft_assembly = "phillippy_plus_cortex_method/sampleDataOxford/windows/polished.fa";
my $draft_assembly = 'phillippy_plus_cortex_method/sampleDataOxford/windows/testassembly.fa';
#my $MISEQ_READS = "/Net/cycloid/data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R1_001_val_1.fq /Net/cycloid/data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R2_001_val_2.fq";
my $base_dir = "phillippy_plus_cortex_method/sampleDataOxford/test/";
my $window_size = 3;
my $NUM_PROCS = 20;
my $k = 31;
my $label = 'test';

# 1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)
print "1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)\n";
#my $cmd = "mkdir -p $base_dir/windows_results";
#my $rcmd = qx{cmd};
#cp phillippy_plus_cortex_method/sampleDataOxford/polished.fa $base_dir
#bwa index $base_dir/polished.fa
#bwa mem $base_dir/polished.fa $MISEQ_READS > $base_dir/miseq_reads.sam
#samtools view -Sb $base_dir/miseq_reads.sam > $base_dir/miseq_reads.bam
#samtools sort -o $base_dir/miseq_reads_sorted.bam $base_dir/miseq_reads.bam
#cp phillippy_plus_cortex_method/sampleDataOxford/windows/""miseq_reads_sorted.bam $base_dir/""miseq_reads_sorted.bam
#samtools index -b $base_dir/""miseq_reads_sorted.bam

# 2. Divide up into windows
print "2. Divide up into $window_size b windows\n";
my $seqio = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	my $contig_read = $seq->seq;
    	my $len_contig = length $contig_read;
	my $par = "parallel --gnu -j $NUM_PROCS \"echo {}; perl -Mpolish_funcs -e \'create_window_and_contents(hi,there);\'\" ::: \$(eval echo {0..$len_contig..$window_size})";
   	#my $par = "parallel --gnu -j $NUM_PROCS \"echo {}; perl -Mpolish_funcs -e \'create_window_and_contents(\'$base_dir\',\"$k\",\"$label\",\"$window_size\");\'\" ::: \$(eval echo {0..$len_contig..$window_size})";
	#my $par = "parallel --gnu -j $NUM_PROCS \"echo {}; perl -Mpolish_funcs -e \'blah\'\" ::: \$(eval echo {0..$len_contig..$window_size})";
	my $rpar = qx{$par};
	print $rpar;
    }

=comment
	contig_read=$(sed -n "/$contig/{n;p;}" $draft_assembly)
	len_contig=$(sed -n "/$contig/{n;p;}" $draft_assembly | wc -c)
	start_pos=0
	while [ $start_pos -lt $len_contig ]; do
		end_pos=$((start_pos + window_size))
		mkdir -p $base_dir/windows_results/$contig.$start_pos\-$end_pos/k$k$label
		contig_read_substring=${contig_read:$start_pos:$window_size}
		echo -e ">$contig.$start_pos-$end_pos\n$contig_read_substring" > $base_dir/windows_results/$contig.$start_pos-$end_pos/k$k$label/draft_assembly.$contig.$start_pos-$end_pos.fa
		samtools view  -b $base_dir/"miseq_reads_sorted.bam" $contig:$start_pos-$end_pos -o $base_dir/windows_results/$contig.$start_pos-$end_pos/k$k$label/miseq_reads.$contig.$start_pos-$end_pos.bam
		samtools bam2fq $base_dir/windows_results/$contig.$start_pos-$end_pos/k$k$label/miseq_reads.$contig.$start_pos-$end_pos.bam > $base_dir/windows_results/$contig.$start_pos-$end_pos/k$k$label/miseq_reads.$contig.$start_pos-$end_pos.fq
		let start_pos=$end_pos
	done
done

# 3. For each window, run_calls and correct
echo "3. For each window, run_calls and correct"
ls -d $PWD/$base_dir/windows_results/*/k$k$label | parallel --gnu -j $NUM_PROCS 'bash code/cortex_correction.sh {} {}/draft_assembly*.fa '$window_size' '$k' {}/*.fq &>>{}/log_cortex_correction_k'$k$label'.txt'

# 4. Combine corrected windows
echo "4. Combine corrected windows"
for contig in $contigs; do
        contig_read=$(sed -n "/$contig/{n;p;}" $draft_assembly)
        len_contig=$(sed -n "/$contig/{n;p;}" $draft_assembly | wc -c)
        start_pos=0
        while [ $start_pos -lt $len_contig ]; do
		if [ $start_pos -eq 0 ]; then
			echo ">$contig" &>> $base_dir/corrected_assembly_k$k$label.fa
		fi
		end_pos=$((start_pos + window_size))
		# want all but first line, and no newlines
		chunk_seq=$(tail -n +2 $base_dir/windows_results/$contig.$start_pos-$end_pos/k$k$label/ref_bc_k$k.fa)
		#len_chunk=`expr length $chunk_seq`
		len_chunk=${#chunk_seq} 
		if [ "$len_chunk" -lt "$window_size" ]; then
			echo "$start_pos seq is $len_chunk bps long, shorter than $window_size"
		fi
		echo "$chunk_seq"  &>> $base_dir/corrected_assembly_k$k$label.fa.tmp
		let start_pos=$end_pos
        done
done
cat $base_dir/corrected_assembly_k$k$label.fa.tmp | tr -d '[:space:]' &>> $base_dir/corrected_assembly_k$k$label.fa
rm $base_dir/corrected_assembly_k$k$label.fa.tmp
=cut
my $footer = "bash code/footer.sh";
my $rfooter = qx{$footer};
