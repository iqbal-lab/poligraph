#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';
use Getopt::Long;

my $poli_dir = abs_path($0);
$poli_dir =~ s/poligraph.pl//;
# parameters needed by this script
my $draft_assembly;
my $reads_fq;
my $base_dir = abs_path('.');
my $window_size = 10000;
my $global;
my $NUM_PROCS = 20;

# parameters passed via polish_window.pl script into cortex_correction.pl script
my $stampy_bin = "/apps/well/stampy/1.0.24-py2.7/stampy.py";
my $vcftools_dir= "vcftools_0.1.13/";
my $cortex_dir = 'cortex/';
my $k = 31;
my $bc = "yes";
my $pd = "no";
my $read_type = "illumina";
my $manual_clean_file = "";
my $auto_clean = "stringent";
my $qthresh = 10;
my $label = "k$k.q$qthresh";
my $genome_size = 5000000;
my $usage = "usage: $0 --draft_assembly path/to/draft_assembly.fa --reads_fq path/to/reads.fq --base_dir dir_to_work_from --window_size int --num_procs int --label string --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py [options]";

my $result = GetOptions (       "draft_assembly=s"              => \$draft_assembly,
				"reads_fq=s"			=> \$reads_fq,
				"base_dir=s"			=> \$base_dir,
                                "window_size=i"                 => \$window_size,
                                "global"			=> \$global,
				"num_procs=i"                   => \$NUM_PROCS,
				"label=s"			=> \$label,
                                "cortex_dir=s"                  => \$cortex_dir,
                                "vcftools_dir=s"                => \$vcftools_dir,
                                "stampy_bin=s"                  => \$stampy_bin,
                                "k=i"                           => \$k,
				"genome_size=s"                 => \$genome_size,
                                "pd=s"                          => \$pd,
                                "bc=s"                          => \$bc,
                                "read_type=s"                   => \$read_type,
                                "manual_clean_file=s"           => \$manual_clean_file,
                                "auto_clean=s"                  => \$auto_clean,
                                "qthresh=i"                     => \$qthresh
                                                ) or die "Incorrect usage. $usage\n";

my $start = "echo \"***** START POLIGRAPH: \$(date)\"";
my $ret_start = qx{$start};
print $ret_start;
my $t1 = "date";
my $t1r = qx{$t1};
print "$t1r\n";

if (!(-d "base_dir"))
{
        my $md = "mkdir -p $base_dir";
        my $ret_md = qx{$md};
}

if ( $global ){
print "\n***** 1. Run cortex correction globally\n";
my $cmd = "perl $poli_dir"."cortex_correction.pl";
$cmd .=" --outdir $base_dir";
$cmd .=" --draft_assembly $draft_assembly";
$cmd .=" --reads \"$reads_fq\"";
$cmd .=" --genome_size $genome_size";
$cmd .=" --cortex_dir $cortex_dir";
$cmd .=" --vcftools_dir $vcftools_dir";
$cmd .=" --stampy_bin $stampy_bin";
$cmd .=" --kmer $k";
$cmd .=" --pd $pd";
$cmd .=" --bc $bc";
$cmd .=" --read_type $read_type";
$cmd .=" --manual_clean_file $manual_clean_file";
$cmd .=" --auto_clean $auto_clean";
$cmd .=" --qthresh $qthresh";
$cmd .=" &>>$base_dir/log_cortex_correction.txt";
print "$cmd\n";
my $rcmd = qx{$cmd};
$t1r = qx{$t1};
print "$t1r\n";

} else {
########################################################################
# 1. Run BWA-MEM to map miseq reads against draft assembly
########################################################################
print "\n***** 1. Run BWA-MEM to map miseq reads against draft assembly\n";

my $cmd1 = "cp $draft_assembly $base_dir/draft_assembly.fa";
print "$cmd1\n";
my $rcmd1 = qx{$cmd1};
my $cmd2 = "bwa index $base_dir/draft_assembly.fa";
print "$cmd2\n";
my $rcmd2 = qx{$cmd2};
my $cmd3 = "bwa mem $base_dir/draft_assembly.fa $reads_fq > $base_dir/reads.sam";
print "$cmd3\n";
my $rcmd3 = qx{$cmd3};
my $cmd4 = "samtools view -Sb $base_dir/reads.sam > $base_dir/reads.bam";
print "$cmd4\n";
my $rcmd4 = qx{$cmd4};
my $cmd5 = "samtools sort $base_dir/reads.bam -o $base_dir/reads_sorted.bam";
print "$cmd5\n";
my $rcmd5 = qx{$cmd5};
my $cmd6 = "samtools index -b $base_dir/reads_sorted.bam";
print "$cmd6\n";
my $rcmd6 = qx{$cmd6};

########################################################################
# 2. Correct in windows
########################################################################
print "\n***** 2. Correct in $window_size b windows\n";
my $seqio = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	my $contig = $seq->id;
    	my $len_contig = length $seq->seq;
	$len_contig = $len_contig - 1;
	my $par = "parallel --gnu -j $NUM_PROCS \"perl $poli_dir"."polish_window.pl";
	$par.=" --contig $contig";
        $par.=" --start_pos {}";
        $par.=" --window_size $window_size";
        $par.=" --base_dir $base_dir";
        $par.=" --draft_assembly $draft_assembly";
        $par.=" --reads_bam $base_dir/reads_sorted.bam";
        $par.=" --stampy_bin $stampy_bin";
        $par.=" --vcftools_dir $vcftools_dir";
        $par.=" --cortex_dir $cortex_dir";
        $par.=" --k $k";
        $par.=" --bc $bc";
	$par.=" --pd $pd";
	$par.=" --read_type $read_type";
	$par.=" --manual_clean_file $manual_clean_file";
	$par.=" --auto_clean $auto_clean";
	$par.=" --qthresh $qthresh";
	$par.=" --label $label";
	$par.="\" ::: \$(eval echo {0..$len_contig..$window_size})";
	#$par.="\" ::: \$(eval echo {0..$window_size..$len_contig})";
	print "$par\n";
	my $rpar = qx{$par};
    }

########################################################################
# 3. Combine corrected windows
########################################################################
print "\n***** 3. Combine corrected windows\n";
#my $seqio2 = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
open (OUTFILE, ">$base_dir/poligraph_corrected.$label.fa");
while(my $seq = $seqio->next_seq) {
	my $contig = $seq->id;
	my $len_contig = length $seq->seq;
	my $corrected_seq = "";
	for (my $start_pos = 0; $start_pos <= $len_contig; $start_pos += $window_size) {
		my $end_pos = $start_pos + $window_size;
		my $seqio2 = Bio::SeqIO->new(-file => "$base_dir/windows/$contig.$start_pos\-$end_pos/$label/poligraph_corrected.fa", '-format' => 'Fasta');
		my $record = $seqio2->next_seq;
		$corrected_seq .= $record->seq;
		print "$end_pos\n";
	}
	print OUTFILE ">$contig\n$corrected_seq";
}
close(OUTFILE);
}
$t1r = qx{$t1};
print "$t1r\n";

my $end = "echo \"***** FINISH POLIGRAPH: \$(date)\"";
my $ret_end = qx{$end};
print $ret_end;
