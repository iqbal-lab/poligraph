#!/user/bin/perl -w
use strict;


my $vcf =  shift;

my $conf_thresh = 100;

my $out = $vcf.".filtered";


open(FILE, $vcf)||die();
open(OUT, ">".$out)||die();

while (<FILE>)
{
    my $line  = $_;
    chomp $line;
    if ($line =~ /^\#/)
    {
	print OUT $line."\n";
    }
    else
    {
	##split out the tab-separated fields into an array
	my @sp = split(/\t/, $line);
	my $filter = $sp[6];
	my $info = $sp[7];
	if ( ($filter eq "PASS") )
	{
		#look at gt_confidence
		my $descriptor = $sp[8]; #   GT:GL:GT_CONF:COV
		##really I should just parse it
		if  ( ($descriptor ne "GT:GL:GT_CONF:COV") && ($descriptor ne "GT:COV:GT_CONF:GL") && ($descriptor ne "GT:COV:GT_CONF") )
		{
		    die("Combined VCF is not in the format I expected, descriptor is $descriptor not GT:GL:GT_CONF:COV GT:COV:GT_CONF:GL or GT:COV:GT_CONF\n");
		}
		my @fields = split(/:/, $sp[9]);
		my $conf=$fields[2];
		my $gt=$fields[0];
		if ( ($conf> $conf_thresh) && ($gt eq "1/1") )
		{
		    print OUT "$line\n";
		}
	}
    }    
}
close(FILE);
close(OUT);
