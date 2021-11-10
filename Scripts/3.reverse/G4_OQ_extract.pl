#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $out, $ca, $lo, $bu, $ot, $tt, $stat, $help);
GetOptions(
	"i=s"=>\$in,
	"o=s"=>\$out,
	"ca=s"=>\$ca,
	"lo=s"=>\$lo,
	"bu=s"=>\$bu,
	"ot=s"=>\$ot,
	"tt=s"=>\$tt,
	"st=s"=>\$stat,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usgae:
	perl $0 [options]
Options:
	-i file,	${red}input oq file$end
	-o file,	${red}output oq file containing G4 structure$end
	-ca file,	${red}output of canonical_unmerged, if want enter the out.name$end
	-lo file,	${red}output of longloop_unmerged, if want enter the out.name$end
	-bu file,	${red}output of bulge_unmerged, if want enter the out.name$end
	-ot file,	${red}output of others_unmerged, if want enter the out.name$end
	-tt file,	${red}output of TwoTrcts_unmerged, if want enter the out.name$end
	-st file,	${red}output of stat results, if want enter the out.name$end
INFO
die $usage if ($help || !$in || !$stat);

if($in =~ /\.gz$/)
{
	open IN,"gzip -dc $in |" || die $!;
}
else
{
	open IN,"$in" || die $!;
}
open OUT,">$out" if $out;
open Ca,">$ca" if $ca;
open TT,">$tt" if $tt;
open Lo,">$lo" if $lo;
open Bu,">$bu" if $bu;
open Ot,">$ot" if $ot;

sub Canonical
{
	my $C = 0;
	$C++ if $_[4] =~ /(C{3}[ATCG]{1,7}?){3}C{3}/;
	return $C;
}

sub Long
{
	my $L = 0;
	my @LongLoop = (
		'C{3}[ATCG]{8,12}?C{3}([ATCG]{1,7}?C{3}){2}',
		'C{3}[ATCG]{1,7}?C{3}[ATCG]{8,21}?C{3}[ATCG]{1,7}?C{3}',
		'(C{3}[ATCG]{1,7}?){2}C{3}[ATCG]{8,12}?C{3}',
		'C{3}[ATCG]{8,12}?C{3}[ATCG]{1,7}?C{3}[ATCG]{8,12}?C{3}'
	);
	for my $i(@LongLoop)
	{
		if($_[4] =~ $i)
		{
			$L ++;
			last;
		}
	}
	return $L;
}

sub Bulge
{
	my $B = 0;
	my @Bulges = (
		'(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)([ATCG]{1,7}?C{3}){3}',
		'C{3}[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)([ATCG]{1,7}?C{3}){2}',
		'(C{3}[ATCG]{1,7}?){2}(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?C{3}',
		'(C{3}[ATCG]{1,7}?){3}(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)',
		'((CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?){2}C{3}[ATCG]{1,7}?C{3}',
		'(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?C{3}[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?C{3}',
		'(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)([ATCG]{1,7}?C{3}){2}[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)',
		'C{3}[ATCG]{1,7}?((CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?){2}C{3}',
		'C{3}[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?C{3}[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)',
		'(C{3}[ATCG]{1,7}?){2}(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)[ATCG]{1,7}?(CC[ATG]{1,7}?C|C[ATG]{1,7}?CC)'
	);
	for my $i(@Bulges)
	{
		if($_[4] =~ $i){
			$B ++;
			last;
		}
	}
	return $B;
}

sub TwoTracts
{
	my $T = 0;
	$T ++ if $_[4] =~ /(C{2}[ATCG]{1,7}?){3}C{2}/;
	return $T;
}

if($stat){
	open ST, "> $stat" || die $!;
}

print ST "Canonical\tLong\tLulge\tTwoTracts\tOthers\n";

my ($canonical, $long, $bulge, $twotracts, $others)=(0, 0, 0, 0, 0);
while(<IN>)
{
	chomp;
	my @temp = split/\s+/, $_;
	if(&Canonical(@temp) > 0)
	{
		print Ca "$_\n" if $ca;
		print OUT "$_\n" if $out;
		$canonical ++;
		next;
	}
	if(&Long(@temp) > 0)
	{
		print Lo "$_\n" if $lo;
		print OUT "$_\n" if $out;
		$long ++;
		next;
	}
	if(&Bulge(@temp) > 0)
	{
		print Bu "$_\n" if $bu;
		print OUT "$_\n" if $out;
		$bulge ++;
		next;
	}
	if(&TwoTracts(@temp) > 0)
	{
		print TT "$_\n" if $tt;
		print OUT "$_\n" if $out;
		$twotracts ++;
		next;
	}
	print Ot "$_\n" if $ot;
	$others ++;
}
close IN;
close OUT;
close Ca if $ca;
close TT if $tt;
close Lo if $lo;
close Bu if $bu;
close Ot if $ot;

print ST "$canonical\t$long\t$bulge\t$twotracts\t$others\n";

close ST if $stat;
