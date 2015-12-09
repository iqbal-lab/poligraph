#!/usr/bin/perl
# Script that identifies reads that map to a region, and runs cortex_correction.pl in that region.

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(min);

# parameters needed by this script
my $contig;
my $start_pos;
my $window_size;
my $base_dir;
my $draft_assembly;
my $reads_bam;

# parameters passed into cortex_correction.pl script
my $stampy_bin = "/apps/well/stampy/1.0.24-py2.7/stampy.py";
my $vcftools_dir= "~/apps/vcftools_0.1.13/";
my $cortex_dir = '~/apps/cortex/';
my $k = 31;
my $bc = "yes";
my $pd = "no";
my $read_type = "illumina";
my $manual_clean_file = "";
my $auto_clean = "stringent";
my $qthresh = 10;
my $label = "k$k.q$qthresh";

my $usage = "usage: $0 --contig contig_id --start_pos int --window_size int --base_dir dir_to_work_from --draft_assembly path/to/draft_assembly.fa --reads_bam path/to/reads.bam --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py [options]";

my $result = GetOptions (       "contig=s"                      => \$contig,
				"start_pos=i"			=> \$start_pos,
				"window_size=i"			=> \$window_size,
				"base_dir=s"			=> \$base_dir,
                                "draft_assembly=s"              => \$draft_assembly,
                                "reads_bam=s"                   => \$reads_bam,
                                "cortex_dir=s"                  => \$cortex_dir,
                                "vcftools_dir=s"                => \$vcftools_dir,
                                "stampy_bin=s"                  => \$stampy_bin,
                                "k=i"				=> \$k,
                                "pd=s"                          => \$pd,
                                "bc=s"                          => \$bc,
                                "read_type=s"                   => \$read_type,
                                "manual_clean_file=s"           => \$manual_clean_file,
                                "auto_clean=s"                  => \$auto_clean,
                                "qthresh=i"                     => \$qthresh,
				"label=s"			=> \$label
                                                ) or die "Incorrect usage. $usage\n";

my $end_pos = $start_pos + $window_size;

my $start = "echo \"***** START WINDOW $contig:$start_pos-$end_pos: \$(date)\n\"";
my $ret_start = qx{$start};
print $ret_start;

# Make dir for this window
my $window_dir = "$base_dir/windows/$contig.$start_pos\-$end_pos/$label";
my $cmd1 = "mkdir -p $window_dir";
my $rcmd1 = qx{$cmd1};
       
# Place draft assembly for window to the dir (i.e. the input to be cleaned)
my $seqio = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	if ( $seq->id eq $contig ){
		my $len_contig = length $seq->seq;
        	my $contig_read_substring = $seq->subseq($start_pos + 1, min($len_contig, $end_pos));
	open (DFILE, ">$window_dir/draft_assembly.$contig.$start_pos-$end_pos.fa");
	print DFILE ">$contig.$start_pos-$end_pos\n$contig_read_substring";
	close(DFILE);
	}
}       

# Make a fastq of miseq reads corresponding to the region
my $cmd2 = "samtools view  -b $reads_bam $contig:$start_pos-$end_pos -o $window_dir/reads.$contig.$start_pos-$end_pos.bam";
#my $rcmd2 = qx{$cmd2};
my $cmd3 = "samtools bam2fq $window_dir/reads.$contig.$start_pos-$end_pos.bam > $window_dir/reads.$contig.$start_pos-$end_pos.fq";
#my $rcmd3 = qx{$cmd3};
        
# Run cortex correction in window
my $cmd4 = "perl cortex_correction.pl";
$cmd4 .=" --outdir $window_dir";
$cmd4 .=" --draft_assembly $draft_assembly";
$cmd4 .=" --reads $window_dir/reads.$contig.$start_pos-$end_pos.fq";
$cmd4 .=" --genome_size $window_size";
$cmd4 .=" --cortex_dir $cortex_dir";
$cmd4 .=" --vcftools_dir $vcftools_dir";
$cmd4 .=" --stampy_bin $stampy_bin";
$cmd4 .=" --kmer $k"; 
$cmd4 .=" --pd $pd"; 
$cmd4 .=" --bc $bc"; 
$cmd4 .=" --read_type $read_type"; 
$cmd4 .=" --manual_clean_file $manual_clean_file"; 
$cmd4 .=" --auto_clean $auto_clean"; 
$cmd4 .=" --qthresh $qthresh"; 
$cmd4 .=" &>>$window_dir/log_cortex_correction.txt";
print "$cmd4\n";
#my $rcmd4 = qx{$cmd4};

my $end = "echo \"***** FINISH WINDOW $contig:$start_pos-$end_pos: \$(date)\n\"";
my $ret_end = qx{$end};
print $ret_end;
