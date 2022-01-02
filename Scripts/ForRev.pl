#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($in, $fai, $outr, $outf, $algo, $help);
GetOptions(
	"i=s"=>\$in,
	"f=s"=>\$fai,
	"of=s"=>\$outf,
	"or=s"=>\$outr,
	"ag=i"=>\$algo,
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
	-ag <INT> <1 represent house Perl method, 0 represent samtools view, 2 represent sambamba view for eatracting reads from inpu file> (dflt=2)
INFO

die $usage if ($help || !$in || !$fai);

die "$red$fai !exists$end\n" if(!(-e "$fai"));

die "$red$in !exists$end\n" if(!(-e "$in"));

if(!(defined $algo)){
	$algo = 2;
}

if($algo == 1){
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
			next unless(!$biFlag[2] && !$biFlag[11]);
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
}
else{
	die "${red}the input must be bam|sam file when use samtools|sambamba view for extract reads$end\n" if($in !~ /\.bam$|\.sam$/);
	if($algo == 2){
		if($in =~ /\.sam$/){
			if($outf){
				my $cmd = "sambamba view -F \"(not reverse_strand) and (not unmapped) and (not supplementary)\" -f bam -l 9 -h -t 20 -S -o $outf $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}

			if($outr){
				my $cmd = "sambamba view -F \"reverse_strand and (not unmapped) and (not supplementary)\" -f bam -l 9 -h -t 20 -S -o $outr $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}
		}
		else{
			if($outf){
				my $cmd = "sambamba view -F \"(not reverse_strand) and (not unmapped) and (not supplementary)\" -f bam -l 9 -h -t 20 -o $outf $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}

			if($outr){
				my $cmd = "sambamba view -F \"reverse_strand and (not unmapped) and (not supplementary)\" -f bam -l 9 -h -t 20 -o $outr $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}
		}
	}
	else{
		if($in =~ /\.sam$/){
			if($outf){
				my $cmd = "samtools view -h -@ 20 --write-index -F 4 -F 2048 -F 16 -o $outf $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}

			if($outr){
				my $cmd = "samtools view -h -@ 20 --write-index -F 4 -F 2048 -f 16 -o $outr $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}
		}
		else{
			if($outf){
				my $cmd = "samtools view -h -@ 20 --write-index -F 4 -F 2048 -F 16 -o $outf $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}

			if($outr){
				my $cmd = "samtools view -h -@ 20 --write-index -F 4 -F 2048 -f 16 -o $outr $in";
				print "\nNOTICE: Running with system command <$cmd>\n";
				system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
			}
		}
	}
}
