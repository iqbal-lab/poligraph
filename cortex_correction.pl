#!/usr/bin/perl

# Script to polish an assembly (or window)  with some reads using cortex
# In stages:
#       1. Make 'ref' graph and stampy hash of the assembly or the window fasta that was input
#       2. Use cortex run_calls to make calls of the reads against the draft assembly (ref graph)
#       3. Update the assembly with corrections that have been called

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

BEGIN
{
    my $cortex_dir = '~/apps/cortex/';
    push ( @INC, $cortex_dir."scripts/calling/");
}

use BasicUtils qw ( check_cortex_runnable add_slash is_fastq count_bases_in_fasta create_dir_if_does_not_exist );

my $outdir;
my $ref_fa;
my $genome_size;
my $k;
my $reads_fq;

my $usage = "usage: $0 --outdir outdir --draft_assembly path/to/draft_assembly.fa --genome_size x --kmer x --reads /path/to/reads.fq";

my $result = GetOptions (       "outdir=s"                   	=> \$sample_name,
                                "draft_assembly=s"              => \$ref_fa,
				"genome_size=i"			=> \$genome_size,
				"kmer=i"			=> \$k,
				"reads=s"			=> \$reads_fq
                                                ) or die "Incorrect usage. $usage\n";
$outdir = BasicUtils::add_slash($outdir);
BasicUtils::create_dir_if_does_not_exist($outdir, "make outdir if doesn't exist");
my $working_dir = abs_path('.');

my $start = "bash code/header2.sh";
my $ret_start = qx{$start};
print $ret_start;

######################################################################################
# 1. Make graph and stampy hash of draft assembly
######################################################################################
print "1. Make graph and stampy hash of draft assembly\n";

BasicUtils::create_dir_if_does_not_exist("$outdir/ref/stampy", "make stampy dir if doesn't exist");
BasicUtils::create_dir_if_does_not_exist("$outdir/ref/ctx_bins", "make cortex dir if doesn't exist");
my $c1 = "ls $ref_fa > $outdir/ref.list";
my $rc1 = qx{$c1};
my $c2 = "cortex_var_31_c1 --kmer_size $k --mem_height 20 --mem_width 100";
$c2 .= " --se_list $outdir/ref.list --max_read_len 10000";
$c2 .= " --dump_binary $outdir/ref/ctx_bins/ref.k$k.ctx --sample_id REF";
my $rc2 = qx{$c2};
my $c3 = "stampy.py -G $outdir/ref/stampy/REF $ref_fa";
my $rc3 = qx{$c3};
my $c4 = "stampy.py -g $outdir/ref/stampy/REF -H $outdir/ref/stampy/REF";
my $rc4 = qx{$c4};

######################################################################################
# 2. Run calls with reads against draft assembly with cortex
######################################################################################
print "2. Run calls with reads against draft assembly with cortex\n";

# Make INDEX for read fastqs
my $c5 = "rm $outdir/reads_fq.list";
my $rc5 = qx{$c5};
my $c6 = "for f in $reads_fq; do echo -e "$(readlink -f $f)" &>> $outdir/reads_fq.list; done";
my $rc6 = qx{$c6};
my $c7 = "echo -e "reads_fq\t$(readlink -f $outdir/reads_fq.list)\t.\t." &> $outdir/INDEX";
my $rc7 = qx{$c7};

# If running in windows, don't need as much memory
if ( $genome_size gt 1000000 ) {mem_height=23}
elsif ( $genome_size lt 5001 ) {mem_height=14}
else {mem_height=16}

# Run cortex run_calls with bubble caller
my $c8 = "perl ~/apps/cortex/scripts/calling/run_calls.pl --first_kmer $k";
$c8 .= " --fastaq_index $outdir/INDEX --manual_override_cleaning $outdir/CLEANFILE";
$c8 .= " --bc yes --pd no --outdir $outdir/results --outvcf bubble_caller_k$k --ploidy 1";
$c8 .= " --stampy_hash $outdir/ref/stampy/REF --stampy_bin /apps/well/stampy/1.0.24-py2.7/stampy.py";
$c8 .= " --list_ref_fasta $outdir/ref.list --refbindir $outdir/ref/ctx_bins/"; 
$c8 .= " --genome_size $genome_size --qthresh 10 --mem_height $mem_height --mem_width 100";
$c8 .= "--vcftools_dir /home/rmnorris/apps/vcftools_0.1.13/ --do_union yes";
$c8 .= " --ref CoordinatesAndInCalling --workflow independent --logfile $outdir/log_run_calls_bc_k$k.txt";
$c8 .= " --gt_assemblies no";
my $rc8 = qx{$c8};

# Run cortex run_calls with path divergence
#my $c9 = "perl ~/apps/cortex/scripts/calling/run_calls.pl --first_kmer $k"; 
#$c9 .= " --fastaq_index $outdir/INDEX --manual_override_cleaning $outdir/CLEANFILE";
#$c9 .= " --bc no --pd yes --outdir $outdir/results --outvcf path_divergence_k$k --ploidy 1";
#$c9 .= " --stampy_hash $outdir/ref/stampy/REF --stampy_bin /apps/well/stampy/1.0.24-py2.7/stampy.py";
#$c9 .= " --list_ref_fasta $outdir/ref.list --refbindir $outdir/ref/ctx_bins/";
#$c9 .= " --genome_size $genome_size --qthresh 10 --mem_height $mem_height --mem_width 100";
#$c9 .= "--vcftools_dir /home/rmnorris/apps/vcftools_0.1.13/ --do_union yes";
#$c9 .= " --ref CoordinatesAndInCalling --workflow independent --logfile $outdir/log_run_calls_pd_k$k.txt";
#$c9 .= " --gt_assemblies no";
#my $rc9 = qx{$c9};

######################################################################################
# 3. Correct the draft assembly with calls
######################################################################################
print "3. Correct the draft assembly with calls\n";

my $c10 = "perl code/filter_cortex_vcf.pl $outdir/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf";
my $rc10 = qx{$c10};
my $c11 = "bgzip $outdir/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered";
my $rc11 = qx{$c11};
my $c12 = "tabix -p vcf $outdir/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz";
my $rc12 = qx{$c12};
my $c13 = "cat $ref_fa | vcf-consensus $outdir/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz > $outdir/ref_bc_k$k.fa 2>$outdir/log_vcf_consensus_k$k.txt";
my $rc13 = qx{$c13};

my $end = "bash code/footer.sh";
my $ret_end = qx{$end};
print $ret_end;
