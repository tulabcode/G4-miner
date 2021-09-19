#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my ($dir, $fdir, @chrs, $help);

GetOptions(
                "d=s"=>\$dir,
                "ch=s"=>\@chrs,
				"f=s"=>\$fdir,
                "help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
        perl $0 [options]
Options:
		-d /mainsd/duanmq/qgrs/GM/OQ/
		-f /mainsd/duanmq/qgrs/GM/PQ/chromosomes
		-ch chr1 -chr chr2
INFO

die $usage if ($help || !$dir || !@chrs);

my @num;
for my $chr (@chrs)
{
	open IN,"$dir/$chr/PQ_g4predict_plus" || die $!;
	open INF,"$fdir/${chr}.fa" || die $!;
	
	my @ref=('N');
	while(<INF>){
		next if $_ =~ />/;
		chomp;
		my @seq=split("",$_);
		push(@ref,@seq);
	}

	while(<IN>){
		chomp;
		my @temp=split/\s+/,$_;
		my $pq_s=$temp[1]-62;
		my $pq_e=$temp[1]+12;
		my $count=0;
		for(my $i=$pq_s;$i<=$pq_e;$i++){
			if(!$ref[$i]){
				last;
			}
			if($ref[$i] eq "G"){
				$count++;
				$num[$count]++;
			}
		}
	}
}

for my$i(0..$#num){
	print "$i\t$num[$i]\n" if $num[$i];
}

close IN;
close INF;
