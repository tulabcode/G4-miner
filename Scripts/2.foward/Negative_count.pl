#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in,$in_median,$out,$region_length,$help);

GetOptions(
		"i=s"=>\$in,
		"m=s"=>\$in_median,	
		"l=i"=>\$region_length,
		"o=s"=>\$out,
		"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file> <input ${red}the negative_sloci.gz or false_positive_sloci.gz$end>
	-m <file> <input ${red}from the output of Median.pl: Median.gz$end>
	-l <int> <the length of detection region>
	-o <file> <output array, ${red}negative_array or false_positive_array$end>
INFO

die $usage if ($help || !$in || !$in_median || !$out || !$region_length);

if($in =~ /\.gz$/)
{
	open IN,"gzip -dc $in |" || die $!;
}
else
{
	open IN,"< $in" || die $!;
}

if($in_median =~ /\.gz$/)
{
	open INM,"gzip -dc $in_median |" || die $!;
}
else
{
	open INM,"< $in_median" || die $!;
}

open OUT,"> $out" || die $!;

my @median = ();
while(<INM>)
{
	chomp;
	my @tempm = split/\s+/, $_;
	$median[$tempm[0]] = $tempm[1];
}
close INM;

my @num;
my $m=12;
my $n=20;
for(my $i = 1; $i <= $n; $i++)
{
	for(my $j = 1; $j <= $m; $j++)
	{
		$num[$i]->[$j] = 0; 
	}
}

sub Max
{
	my $max=0;
	foreach my $element(@_)
	{
		$max = ($element>$max) ? $element : $max;
	}
	return $max;
}

while(<IN>)
{
	chomp;
	my @temp = split/\s+/, $_;
	my @sub_seq = @median[$temp[0]..($temp[0] + $region_length - 1)];
	my @flag_array = @median[($temp[0] + $region_length)..($temp[0] + $region_length + 9)];
	my $flag = &Max(@flag_array);
	for(my $step = 1; $step <= 12; $step++)
	{
		my $count = 0;
		foreach my $loci(@sub_seq)
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

for(my $i = 1; $i <= $n; $i++)
{
	for(my $j = 1; $j <= $m; $j++)
	{
		print OUT "$num[$i]->[$j]\t";
	}
	print OUT "\n";
}
close OUT;
