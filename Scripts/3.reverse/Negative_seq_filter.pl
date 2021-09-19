#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $out, $fasta, $c, $region_length, $help);
GetOptions(
	"i=s"=>\$in,
	"f=s"=>\$fasta,
	"o=s"=>\$out,
	"l"=>\$region_length,
	"c=f"=>\$c,
	"help|?"=>\$help,
);
my $usage=<<INFO;
Usgae:
	perl $0 [options]
Options:
	-i file	<${red}input negative loci gz file of a special chromosome$end>
	-f file	<${red}input fasta file of special chromosome$end>
	-o file	<${red}output false negative loci file$end>
	-l int <${red}length of detection region$end>
	-c float <${red}cutoff of C ratio$end>
INFO
die $usage if ($help || !$in || !$out || !$fasta || !$c || !$region_length);

if($in =~ /\.gz$/)
{
	open IN, "gzip -dc $in |" || die $!;
}
else
{
	open IN, "$in" || die $!;
}

open IN2,"$fasta" || die $!;

if($out =~ /\.gz$/)
{
	open OUT, "| gzip > $out" || die $!;
}
else
{
	open OUT, "| gzip > $out.gz" || die $!;
}

my @seq = ();
while(<IN2>)
{
	chomp;
	next if $_ =~ />/;
	my @line = split("", $_);
	push(@seq, @line);
}
close IN2;

while(<IN>)
{
	chomp;
	my $line = join("", @seq[$_..($_ + $region_length - 1)]);
	my $count = 0;
	$count ++ while($line=~/C/g);
	print OUT "$_\n" if $count>=($region_length*$c);
}
close IN;
close OUT;
