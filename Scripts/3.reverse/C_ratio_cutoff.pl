#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in,$out,$c,$help);
GetOptions(
	"i=s"=>\$in,
	"o=s"=>\$out,
	"c=f"=>\$c
);
my $usage=<<INFO;
Usgae:
	perl $0 [options]
Options:
	-i file,  ${red}input oq file$end
	-o file,  ${red}output OQ_unmerged$end
	-c float, ${red}the cutoff of C ratio, different specie with different cutoff$end
INFO
die $usage if ($help || !$in || !$out || !$c);

if($in =~ /\.gz$/)
{
	open IN,"gzip -dc $in |" || die $!;
}
else
{
	open IN,"< $in" || die $!;
}
open OUT," > $out" ||die $!;

while(<IN>)
{
	chomp;
	my @temp = split/\s+/, $_;
	my $count = 0;
	$count++ while($temp[3]=~/C/g);
	print OUT "$_\n" if $count >= (length($temp[3])*$c);
}
close IN;
close OUT;
