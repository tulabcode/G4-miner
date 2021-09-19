#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in,$dir,$ref,$chrs,$help);

GetOptions(
	"i=s"=>\$in,
	"d=s"=>\$dir,
	"cs=s"=>\$chrs,
	"g=s"=>\$ref,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usgae:
        perl $0 [options]
Options:
        -i <file> <input oq file, ${red}the name of oq file${end}>
        -d <file> <the directory of each chromosome oq file, for eaxmple: ${red}/mainsd/duanmq/qgrs/GM/OQ/${end}chr*/\$in>
        -cs <file> <${red}the all chromosomes of this sample, for example: "chr1,chr2,chr...*,chrX"$end>
		-g <file> <${red}the fasta file of the whole genome reference$end>
		$red-g | -cs must have one exists$end
INFO
die $usage if ($help || !$in || !$dir || !$chrs);

die $usage if (!$ref && !$chrs);

open OUT,"> $dir/$in" || die $!;

my @num;
my $m = 12;
my $n = 20;
for(my $i = 1; $i <= $n; $i++)
{
	for(my $j = 1; $j <= $m; $j++)
	{
		$num[$i]->[$j] = 0;
	}
}

my @chrs;
if($chrs)
{
	@chrs = ();
	$chrs=~s/^\s+|\s+$//g;
	@chrs = split /,/, $chrs;
}
else
{
	@chrs = ();
	die $red,"$ref !exists\n$end" if(!(-e "$ref"));
	open IF,"grep '>' $ref |" || die $!;
	while(<IF>)
	{
		chomp;
		my ($i) = $_ =~ />(.*)/;
		push @chrs, $i;
	}
	close IF;
}


for my $i(@chrs)
{
	my $file;
	$file="$dir/$i/$in";
	open IN,"< $file" || die "cannot open this file!\n";
	my $row = 1;
	while(<IN>)
	{
		chomp;
		my @temp = split/\s+/, $_;
		for my $j(0..$#temp)
		{
			$num[$row]->[$j+1] += $temp[$j];
		}
		$row++;
	}
	close IN;
}

for(my $i = 1; $i <= $n; $i++)
{
	for(my $j = 1; $j <= $m; $j++)
	{
		print OUT "$num[$i]->[$j]\t";
	}
	print OUT "\n";
}
close OUT;
