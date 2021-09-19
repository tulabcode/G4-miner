#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $in_median, $out, $region_length, $help);

GetOptions(
		"i=s"=>\$in,
		"m=s"=>\$in_median,
		"o=s"=>\$out,
		"l=i"=>\$region_length,
		"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file> <input PQs>
	-m <file> <input Median.gz>
	-l <int> <the length of detection region, ${red}defalut=75nt$end>
	-o <file> <output file for saving positive array>
INFO
die $usage if ($help || !$in || !$in_median || !$out);

if($in_median =~ /\.gz$/)
{
	open INM, "gzip -dc $in_median |" || die $!;
}
else
{
	open INM, "$in_median" || die $!;
}

open IN, "$in" || die $!;
open OUT,"> $out" || die $!;

$region_length ||= 75;

my @median = ();
while(<INM>)
{
	chomp;
	my @tempm = split/\s+/, $_;
	$median[$tempm[0]] = $tempm[1];
}
close INM;

my @num = ();
my $m = 12;
my $n = 20;
for(my $i = 1; $i <= $n; $i ++)
{
	for(my $j = 1; $j <= $m; $j ++)
	{
		$num[$i]->[$j]=0;
	}
}

sub Max
{
	my $max = 0;
	foreach my $element(@_)
	{
		$max = ($element > $max) ? $element : $max;
	}
	return $max;
}

while(<IN>)
{
	chomp;
	my @temp = split/\s+/,$_;
	next if(!$median[$temp[0]] or !$median[$temp[1]]);
	my @sub_pq = ();
	my $flag_s = $temp[1] - 12;
	for(my $i = ($flag_s + $region_length); $i > $flag_s; $i --)
	{
		if($median[$i])
		{
			push(@sub_pq, $median[$i]);
		}
		else{
			push(@sub_pq, 0);
		}
	}
	my @flag_array = ();
	for(my $i = 0; $i <= 9; $i ++)
	{
		if($median[$flag_s - $i])
		{
			push(@flag_array, $median[$flag_s - $i]);
		}
		else
		{
			push(@flag_array, 0);
		}
	}
	my $flag = &Max(@flag_array);
	for(my $step = 1; $step <= 12; $step ++)
	{
		my $count = 0;
		foreach my $loci(@sub_pq)
		{
			my $quality = $flag - $step;
			if($loci <= $quality)
			{
				$count ++;
				$num[$count]->[$step] ++;
			}
		}
	}
}
close IN;

for(my $i = 1; $i <= $n; $i ++)
{
	for(my $j = 1; $j <= $m; $j ++)
	{
		print OUT "$num[$i]->[$j]\t";
	}
	print OUT "\n";
}
close OUT;
