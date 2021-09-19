#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $fai, $outr, $outf, $help);
GetOptions(
	"i=s"=>\$in,
	"f=s"=>\$fai,
	"of=s"=>\$outf,
	"or=s"=>\$outr,
	"help|?"=>\$help
);
my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
	-i <file> <${red}input .bam/.sam mapped file$end>
	-f <file> <${red}the fai of the whole genome fasta$end>
	-or <file> <${red}output reverse reads .bam file if demand$end>
	-of <file> <${red}output forward reads .bam file if demand$end>
INFO
die $usage if ($help || !$in || !$fai);

die "$red$fai !exists$end\n" if(!(-e "$fai"));

die "$red$in !exists$end\n" if(!(-e "$in"));

if($in =~ /\.bam$/)
{
	open IN,"samtools view -h $in |" || die $!;
}
elsif($in =~ /\.gz$/)
{
	open IN,"gzip -dc $in |" || die $!;
}
else
{
	open IN,"$in" || die $!;
}

if($outf)
{
	open OUTF, "| samtools view -Sb -@ 10 -t $fai -o $outf - " || die $!;
}
if($outr)
{
	open OUTR, "| samtools view -Sb -@ 10 -t $fai -o $outr - " || die $!;
}

while(<IN>)
{
	chomp;
	if(/^\@/)
	{
		print OUTR "$_\n" if $outr;
		print OUTF "$_\n" if $outf;
	}
	else
	{
		my @temp = split/\t/, $_;
		my $flag = reverse sprintf("%b", $temp[1]);
		my @biFlag = split //,$flag;
		if($biFlag[4])
		{
			print OUTR "$_\n" if $outr;
		}
		else
		{
			print OUTF "$_\n" if $outf;
		}
	}
}
close IN;
close OUTF if $outf;
close OUTR if $outr;

if($outf)
{
	my $cmd = "samtools index -@ 20 $outf";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
}

if($outr)
{
	my $cmd = "samtools index -@ 20 $outr";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
}
