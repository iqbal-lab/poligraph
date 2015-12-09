#!/usr/bin/perl
use strict;
use warnings;

#sub create_window_and_contents ($base_dir, $k, $label, $window_size, $seq, $start_pos) {
sub create_window_and_contents {
	my ($base_dir, $k, $label, $window_size, $seq, $start_pos) = @_;
	print @_;
	#print "$base_dir\n";
	print "$start_pos\n";
#        my $end_pos = $start_pos + $WINDOW_SIZE;
#        my $contig = $seq->id;
#
#        # Make dir for this window
#        my $window_dir = "$DIR/windows_results/$contig.$start_pos\-$end_pos/k$k$label";
#        my $cmd1 = "mkdir -p $window_dir"
#        my $rcmd1 = qx{$cmd1};
#       
#        # Place draft assembly for window to the dir (i.e. the input to be cleaned)
#        my $contig_read_substring = $seq->subseq($start_pos + 1,$end_pos + 1);
#       	open (DFILE, ">>$window_dir/draft_assembly.$contig.$start_pos-$end_pos.fa");
#        print DFILE ">$contig.$start_pos-$end_pos\n$contig_read_substring";
#        close(DFILE);
#       
#        # Make a fastq of miseq reads corresponding to the region
#        my $cmd2 = "samtools view  -b $DIR/miseq_reads_sorted.bam $contig:$start_pos-$end_pos -o $window_dir/miseq_reads.$contig.$start_pos-$end_pos.bam"
#        my $rcmd2 = qx{$cmd2};
#        my $cmd3 = "samtools bam2fq $windows_dir/miseq_reads.$contig.$start_pos-$end_pos.bam > $windows_dir/miseq_reads.$contig.$start_pos-$end_pos.fq"
#        my $rcmd3 = qx{$cmd3};
#        
#        # Run cortex correction in window
#        my $cmd4 = "bash code/cortex_correction.sh $window_dir $window_dir/draft_assembly*.fa $WINDOW_SIZE $k $window_dir/*.fq &>>$window_dir/log_cortex_correction.txt";
}

sub blah {
    print "Ahhh\n";
}
return 1
