#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $in_median, $region_length, $out, $help);

GetOptions(
		"i=s"=>\$in,
		"m=s"=>\$in_median,
		"l=i"=>\$region_length,
		"o=s"=>\$out,
		"help|?"=>\$help
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file> <${red}input PQs_g4predict$end>
	-m <file> <${red}input Median.gz$end>
	-l <int> <the length of detection region$end>
	-o <file> <${red}output start loci of negative reads, .gz$end>
INFO
die $usage if ($help || !$in || !$in_median || !$out);

if($in_median =~ /\.gz$/)
{
	open INM, "gzip -dc $in_median |" || die $!;
}
else
{
	open INM, "< $in_median" || die $!;
}

open IN, "$in" || die $!;

if($out =~ /\.gz$/)
{
	open OUT, "| gzip > $out " || die $!;
}
else
{
	open OUT, "| gzip > $out.gz " || die $!;
}

my @median = ();
while(<INM>)
{
	chomp;
	my @tempm = split/\s+/, $_;
	$median[$tempm[0]] = $tempm[1];
}
close INM;

my $current = 0;
my $pq_e = 0;
while(<IN>)
{
	chomp;
	my @temp = split/\s+/,$_;
	for($current = ($pq_e + 1); $current <= ($temp[0] - 150 - $region_length); $current += $region_length)
	{
		foreach my $i($current..$#median)
		{
			if($median[$i])
			{
				$current = $i;
				last;
			}
		}
		my $consist = 0;
		for(my $i = -10; $i <= ($region_length - 1); $i ++)
		{
			my $loci = $current + $i;
			if(!$median[$loci])
			{
				$consist = 0;
				last;
			}
			$consist ++;
		}
		if($consist)
		{
			print OUT "$current\n";
		}
	}
	$pq_e = $temp[1] + 150;
}
close IN;
close OUT;
