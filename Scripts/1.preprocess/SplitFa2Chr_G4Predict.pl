#!/usr/bin/perl -w
#########################################################################
# FileName: SplitFa2Chr_G4Predict.pl
# Version: 96a78e56-2c37-4de2-9a77-1ce128c2d58e
# Author: Luna <nlu@seu.edu.cn>
# CreatedTime: Thu Sep  9 16:58:20 2021
#########################################################################
use strict;
use Getopt::Long;
use File::Basename;
my $red = "\033[0;31m";
my $end = "\033[0m";
my ($in, $dir, $py2, $G4Predict, $help);

GetOptions
(
	"i|in=s"=>\$in,
	"d|dir=s"=>\$dir,
	"p|py=s"=>\$py2,
	"g|g4=s"=>\$G4Predict,
	"help|?"=>\$help,
);


my $usage=<<INFO;
Usage:
	${red}1. Split the whole genome reference fasta to each chr fasta$end
	${red}2. Using g4predict intra for predicting putative G Quadruplexes in each chromosome$end

	perl $0 
		-i <the whole genome reference fasta file> 
		${red}-d <the directory for saving the chromosome fa file> $end
		-p <the absolute path of python2> (default: `which python2`)
		-g <the absolute path of g4predict> (default: `which g4predict`)
INFO

die "$usage\n" if ($help || !$in || !$dir);

$py2 ||= `which python2`;
$py2=~s/^\s+|\s+$//g;
die "${red}There did not install python2 $py2$end\n" if(!$py2);

$G4Predict ||= `which g4predict`;
$G4Predict=~s/^\s+|\s+$//g;
die "${red}There did not install g4predict $G4Predict$end\n" if(!$G4Predict);

if($in =~ /\.gz$/)
{
	open IN, "gzip -dc $in |" || die $!;
}
else
{
	open IN,"$in" || die $!;
}

my @chrs = ();

while (<IN>)
{
	chomp;
	if ($_ =~ />(.*)/)
	{
		my $chr = $1;

		my $out = "$chr.fa";
		if(!-e "$dir/$chr")
		{
			`mkdir -p $dir/$chr` || die "${red}Error: can't create $dir/$chr$end\n";
		}
		push(@chrs, $chr);
		open OUT, "> $dir/$chr/$out" || die $!;
		print OUT "$_\n";
	}
	else
	{
		print OUT "$_\n";
	}
}
close IN;
close OUT;

for my $chr(@chrs)
{
	my $fa = "$dir/$chr/$chr.fa";
	my $cmd = "$py2 $G4Predict intra -f $fa -b $dir/$chr/PQ_g4predict.bed -s -M";

	print "NOTICE: run this command\n${red}$cmd$end\n";
	system('$cmd');
}

print "${red}$in was splitted into chromosome fa in $dir$end\n";
print "${red}PQs of each chromosome were predicted by $G4Predict, saved into $dir/chr*$end\n";
