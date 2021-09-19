#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $fasta, $out, $cut, $qn, $help);

GetOptions(
		"i=s"=>\$in,
		"f=s"=>\$fasta,
		"o=s"=>\$out,
		"c=i"=>\$cut,
		"q=i"=>\$qn,
		"help|?"=>\$help
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file> <${red}input Median.gz file of special chromosome$end>
	-f <file> <${red}fasta file of special chromosome$end>
	-o <file> <${red}output G4 positions and sequences (commpressed file with .gz)$end>
	-c <int> <${red}the difference between high and low Phred score$end>
	-q <int> <${red}the number that should meet the threshold above, default = 0$end>
INFO
die $usage if ($help || !$in || !$fasta || !$out || !$cut);

if($in =~ /\.gz$/)
{
	open IN, "gzip -dc $in |" || die $!;
}
else
{
	open IN, "$in" || die $!;
}

open INF, "$fasta" || die $!;

if($out =~ /\.gz$/)
{
	open OUT,"| gzip > $out" || die $!;
}
else
{
	open OUT,"| gzip > $out.gz" || die $!;
}

my @ref = ();
while(<INF>)
{
	next if $_ =~ />/;
	chomp;
	my @line = split("", $_);
	push(@ref, @line);
}
close INF;

my @seq = ();
while(<IN>)
{
	chomp;
	my @temp = split/\s+/, $_;
	$seq[$temp[0]] = $temp[1];
}
close IN;

my $count = 0;
$qn ||= 0;

for(my $i = 0; $i <= $#seq; $i ++)
{
	if(!$seq[$i] or !$seq[$i - 1])
	{
		next;
	}
	if(($seq[$i - 1] - $cut) >= $seq[$i])
	{
		my $qualified_num = 1;
		for(my $j = $i + 0; $j <= $i + 74; $j ++)
		{
			if(!$seq[$j])
			{
				$qualified_num = 0;
				last;
			}
			else
			{
				$qualified_num ++ if(($seq[$i - 1] - $cut) >= $seq[$j]);
			}
		}
		if($qualified_num >= $qn)
		{
			my $oqs = $i - 35;
			my $oqe = $i + 74;
			$count ++;
			if($ref[$oqs] and $ref[$oqe])
			{
				my $oq = join('', @ref[$oqs..$oqe]);
				print OUT "$count\t$oqs\t$oqe\t$oq\n";
			}
		}
	}
}
close OUT;
