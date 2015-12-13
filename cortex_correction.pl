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

#### Dear User - you may want to edit stampy_bin, vcftools_dir and cortex_dir
my $stampy_bin = "/apps/well/stampy/1.0.24-py2.7/stampy.py";
my $vcftools_dir= "~/apps/vcftools_0.1.13/";
my $cortex_dir = '~/apps/cortex/';

my $outdir;
my $ref_fa;
my $genome_size;
my $k = 31;
my $reads_fq;
my $bc = "yes";
my $pd = "no";
my $read_type = "illumina";
my $manual_clean_file = "";
my $auto_clean = "stringent";
my $qthresh = 10;

my $usage = "usage: $0 --outdir outdir --draft_assembly path/to/draft_assembly.fa --reads path/to/reads.fq --genome_size genome_size --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py [options]";

my $result = GetOptions (       "outdir=s"                   	=> \$outdir,
                                "draft_assembly=s"              => \$ref_fa,
				"reads=s"			=> \$reads_fq,
				"genome_size=i"                 => \$genome_size,
				"cortex_dir=s"			=> \$cortex_dir,
				"vcftools_dir=s"		=> \$vcftools_dir,
				"stampy_bin=s"			=> \$stampy_bin,
				"kmer=i"                        => \$k,
				"pd=s"				=> \$pd,
				"bc=s"				=> \$bc,
				"read_type=s"			=> \$read_type,
				"manual_clean_file=s"		=> \$manual_clean_file,
				"auto_clean=s"			=> \$auto_clean,
				"qthresh=i"			=> \$qthresh
                                                ) or die "Incorrect usage. $usage\n";
my $working_dir = abs_path('.');

my $start = "echo \"***** START CORRECTION: \$(date)\n\"";
my $ret_start = qx{$start};
print $ret_start;

######################################################################################
# 1. Make graph and stampy hash of draft assembly
######################################################################################
print "***** 1. Make graph and stampy hash of draft assembly\n";
if (!(-d "$outdir/ref/stampy"))
{
        my $md = "mkdir -p $outdir/ref/stampy";
        my $ret_md = qx{$md};
}
if (!(-d "$outdir/ref/ctx_bins"))
{
        my $md2 = "mkdir -p $outdir/ref/ctx_bins";
        my $ret_md2 = qx{$md2};
}

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
print "***** 2. Run calls with reads against draft assembly with cortex\n";

# Make INDEX for read fastqs
if ( -e "$outdir/reads_fq.list" )
{
	my $c5 = "rm $outdir/reads_fq.list";
	my $rc5 = qx{$c5};
}
my $c6 = "for f in $reads_fq; do echo -e \$(readlink -f \$f) &>> $outdir/reads_fq.list; done";
print $c6;
my $rc6 = qx{$c6};
my $c7 = "echo -e \"reads_fq\t\$(readlink -f $outdir/reads_fq.list)\t.\t.\" &> $outdir/INDEX";
my $rc7 = qx{$c7};

# If running in windows, don't need as much memory
my $mem_height;
if ( $genome_size gt 1000000 ) {$mem_height=23;}
elsif ( $genome_size lt 5001 ) {$mem_height=14;}
else {$mem_height=16;}

# Put together the cleaning command from the command line options input
my $clean_cmd;
if ( $manual_clean_file ne "" ) {
	$clean_cmd = "--manual_override_cleaning $manual_clean_file";
} else {
	$clean_cmd = "--auto_clean $auto_clean";
}

# Run cortex run_calls
my $c8 = "perl $cortex_dir/scripts/calling/run_calls.pl --first_kmer $k";
$c8 .= " --fastaq_index $outdir/INDEX $clean_cmd";
$c8 .= " --bc $bc --pd $pd --outdir $outdir/results --outvcf output_k$k --ploidy 1";
$c8 .= " --stampy_hash $outdir/ref/stampy/REF --stampy_bin $stampy_bin";
$c8 .= " --list_ref_fasta $outdir/ref.list --refbindir $outdir/ref/ctx_bins/"; 
$c8 .= " --genome_size $genome_size --qthresh $qthresh --mem_height $mem_height --mem_width 100";
$c8 .= " --vcftools_dir $vcftools_dir --do_union yes";
$c8 .= " --ref CoordinatesAndInCalling --workflow independent --logfile $outdir/log_run_calls_k$k.txt";
$c8 .= " --gt_assemblies no";
print "$c8\n";
my $rc8 = qx{$c8};

######################################################################################
# 3. Correct the draft assembly with calls
######################################################################################
print "***** 3. Correct the draft assembly with calls\n";

my $c9 = "perl filter_cortex_vcf.pl $outdir/results/vcfs/output_k".$k."_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf";
my $rc9 = qx{$c9};
my $c10 = "bgzip $outdir/results/vcfs/output_k".$k."_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered";
my $rc10 = qx{$c10};
my $c11 = "tabix -p vcf $outdir/results/vcfs/output_k".$k."_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz";
my $rc11 = qx{$c11};
my $c12 = "cat $ref_fa | vcf-consensus $outdir/results/vcfs/output_k".$k."_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz > $outdir/poligraph_corrected.fa 2>$outdir/log_vcf_consensus_k$k.txt";
my $rc12 = qx{$c12};

my $end = "echo \"***** FINISH CORRECTION: \$(date)\"";
my $ret_end = qx{$end};
print $ret_end;
