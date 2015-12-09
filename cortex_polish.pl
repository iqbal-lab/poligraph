#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';

# parameters needed by this script
my $draft_assembly;
my $reads_fq;
my $base_dir = abs_path('.');
my $window_size = 10000;
my $NUM_PROCS = 20;

# parameters passed via polish_window.pl script into cortex_correction.pl script
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

my $usage = "usage: $0 --draft_assembly path/to/draft_assembly.fa --reads_fq path/to/reads.fq --base_dir dir_to_work_from --window_size int --num_procs int --label string --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py [options]";

my $result = GetOptions (       "draft_assembly=s"              => \$draft_assembly,
				"reads_fq=s"			=> \$reads_fq,
				"base_dir=s"			=> \$base_dir,
                                "window_size=i"                 => \$window_size,
                                "num_procs=i"                   => \$NUM_PROCS,
				"label=s"			=> \$label,
                                "cortex_dir=s"                  => \$cortex_dir,
                                "vcftools_dir=s"                => \$vcftools_dir,
                                "stampy_bin=s"                  => \$stampy_bin,
                                "k=i"                           => \$k,
                                "pd=s"                          => \$pd,
                                "bc=s"                          => \$bc,
                                "read_type=s"                   => \$read_type,
                                "manual_clean_file=s"           => \$manual_clean_file,
                                "auto_clean=s"                  => \$auto_clean,
                                "qthresh=i"                     => \$qthresh
                                                ) or die "Incorrect usage. $usage\n";


########################################################################
# 1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)
########################################################################
print "1. Run BWA-MEM to map miseq reads against draft assembly (nanopolished)\n";
my $cmd1 = "cp $draft_assembly $base_dir/draft_assembly.fa";
my $cmd2 = "bwa index $base_dir/draft_assembly.fa";
my $cmd3 = "bwa mem $base_dir/draft_assembly.fa $reads_fq > $base_dir/reads.sam";
my $cmd4 = "samtools view -Sb $base_dir/reads.sam > $base_dir/reads.bam";
my $cmd5 = "samtools sort -o $base_dir/reads_sorted.bam $base_dir/reads.bam";
my $cmd6 = "samtools index -b $base_dir/reads_sorted.bam";

########################################################################
# 2. Correct in windows
########################################################################
print "2. Correct in $window_size b windows\n";
my $seqio = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	my $contig = $seq->id;
    	my $len_contig = length $seq->seq;
	my $par = "parallel --gnu -j $NUM_PROCS \"echo {}; perl polish_window.pl";
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
	print "$par\n";
	my $rpar = qx{$par};
	print $rpar;
    }

########################################################################
# 3. Combine corrected windows
########################################################################
print "2. Combine corrected windows\n";
my $seqio2 = Bio::SeqIO->new(-file => $draft_assembly, '-format' => 'Fasta');
open (OUTFILE, ">$base_dir/poligraph_corrected.$label.fa");
while(my $seq = $seqio2->next_seq) {
	my $contig = $seq->id;
	my $len_contig = length $seq->seq;
	my $corrected_seq = "";
	for (my $start_pos = 0; $start_pos <= $len_contig; $start_pos += $window_size) {
		my $end_pos = $start_pos + $window_size;
		my $seqio3 = Bio::SeqIO->new(-file => "$base_dir/windows/$contig.$start_pos\-$end_pos/$label/poligraph_corrected.fa", '-format' => 'Fasta');
		my $record = $seqio3->next_seq;
		$corrected_seq .= $record->seq;
	}
	print OUTFILE ">$contig\n$corrected_seq";
}
close(OUTFILE);
