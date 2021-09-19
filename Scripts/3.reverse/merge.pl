#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $out, $fasta, $help);
GetOptions(
	"i=s"=>\$in,
	"f=s"=>\$fasta,
	"o=s"=>\$out,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usgae:
        perl $0 [options]
Options:
        -i file <${red}input oq file$end>
        -f file <${red}input fasta file$end>
        -o file output file
INFO
die $usage if ($help || !$in || !$out || !$fasta);

if($in =~ /\.gz$/)
{
	open IN, "gzip -dc $in |" || die $!;
}
else
{
	open IN, "$in" || die $!;
}

open IN2, "$fasta" || die $!;
open OUT, ">$out" || die $!;

my @seq = ();
while(<IN2>)
{
	chomp;
	next if $_ =~ />/;
	my @line = split("", $_);
	push(@seq, @line);
}
close IN2;

my $oqs = 0;
my $oqe = 0;
my $num = 0;
while(<IN>)
{
	next if $_ =~ /N/;
	chomp;
	my @temp = split/\s+/, $_;
	if($temp[1] > $oqe)
	{
		if($oqs != 0)
		{
			my $oq = join('', @seq[$oqs..$oqe]);
			$num ++;
			my $oqlen = $oqe - $oqs + 1;
			print OUT "$num\t$oqs\t$oqe\t$oqlen\t$oq\n";
		}
		$oqs = $temp[1];
	}
	$oqe = ($oqe > $temp[2]) ? $oqe : $temp[2];
	$oqe = ($oqe < $#seq) ? $oqe : $#seq;
}
close IN;
close OUT;
