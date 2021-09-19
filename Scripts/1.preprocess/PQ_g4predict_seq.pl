#!usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use File::Basename;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in,$out,$fasta,$help);

GetOptions(
	"i=s"=>\$in,
	"f=s"=>\$fasta,
	"o=s"=>\$out
);

my $usage=<<INFO;
Usgae:
	perl $0 [options]
Options:
	-i file, input PQ file, the directory of $red\$in will be used for save the output$end
	-f file, input fasta file
	-o file, the prefix of the output PQ with seq file, ${red}default [PQ_g4predict] [the output is \$out_plus && \$out_minus]$end
INFO

die $usage if ($help || !$in || !$fasta);

my $dir = dirname(File::Spec->rel2abs( $in ));
$dir=~s/^\s+|\s+$//g;

$out ||= "PQ_g4predict";

open IN,"<$in" || die $!;
open OUTP,">$dir/${out}_plus" || die $!;
open OUTM,">$dir/${out}_minus" || die $!;
open INF,"<$fasta" || die $!;

my @ref=('N');
while(<INF>)
{
	next if $_ =~ />/;
	chomp;
	my @fasta=split("",$_);
	push(@ref,@fasta);
}
close INF;

while(<IN>)
{
	chomp;
	my @temp = split/\s+/,$_;
	my $pqs = $temp[1] + 1;
	my $pqe = $temp[2];
	my $line = join('',@ref[$pqs..$pqe]);
	my $length = length($line);
	if($temp[5] eq "+")
	{
		print OUTP "$pqs\t$pqe\t$length\t$line\t$temp[4]\t$temp[5]\n";
	}
	elsif($temp[5] eq "-")
	{
		print OUTM "$pqs\t$pqe\t$length\t$line\t$temp[4]\t$temp[5]\n";
	}
}
close IN;
close OUTP;
close OUTM;
