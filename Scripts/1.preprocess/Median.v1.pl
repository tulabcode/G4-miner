#!usr/bin/perl -w
use strict;
use Getopt::Long;

my $red = "\033[0;31m";
my $end = "\033[0m";

my($in,$out,$len,$chr,$help);

GetOptions
(
	"i=s"=>\$in,
	"o=s"=>\$out,
	"l=i"=>\$len,
	"ch=s"=>\$chr,
	"help|?"=>\$help
);

my $usage=<<INFO;
Notice:
	${red}The input .sam/.bam file must be sorted according to coordinate.$end
Usage:
	perl $0	[options]
Options:
	-i <string> <in.input sam/bam file>
	-o <string> <out.output median of each positions>
	-l <int> <the longest length of reads, default=150>
	-ch <string> <the number of chromosome>
INFO

die $usage if ($help || !$in || !$out ||!$chr);

die "$red$in does not exists$end\n" if(!(-e "$in"));

if($in =~ /\.bam$/)
{
	my $head = `samtools view -H $in | grep SO`;
	$head =~ s/^\s+|\s+$//g;
	die "$red$in is not sorted bam file$end\n" if($head !~ /coordinate/);
	my $bai=$in;
	$bai=~s/\.bam/\.bai/g;
	if( (!-e "${in}.bai") || (!-e "$bai") )
	{
		`samtools index -@ 20 $in`;
	}
	open IN,"samtools view $in $chr |" || die $!;
}
else
{
	my $so = &Sorted($in);
	die "$red$in is not sorted sam file$end\n" if($so == 0);
	if($in =~ /\.gz$/)
	{
		open IN,"gzip -dc $in | awk '$3==$chr' | " || die $!;
	}
	else
	{
		open IN,"grep -w $chr $in |" || die $!;
	}
}

sub Sorted
{
	my $in = @_;
	if($in =~ /\.gz/)
	{
		open IN, "gzip -dc $in | grep -v '^\@' | " || die $!;
	}
	else
	{
		open IN, "$in | grep -v '^\@' | " || die $!;
	}
	
	my $read = 0;
	my $pos = 0;
	my $tags = 0;
	my $tagu = 0;
	my $chrt;
	while(my $line = <IN>)
	{
		chomp $line;
		next unless($line !~ /^@/);
		$read ++;
		last if($read > 2000);
		my @F = split /\t/, $line;
		if($read == 1)
		{
			$pos = $F[3];
			$chrt = $F[2]; 
		}
		else
		{
			last if($F[2] ne $chrt);
 			if($F[3] >= $pos)
			{
				$tags ++;
				$pos = $F[3];
				$chrt = $F[2];
			}
			else
			{
				$tagu ++;
				$pos = $F[3];
				$chrt = $F[2];
			}
		}
	}
	
	my $sort = ($tags > $tagu && $tagu == 0) ? 1 : 0;
	return $sort;
}

if($out =~ /\.gz$/)
{
	open OUT,"| gzip > $out" || die $!;
}
else
{
	open OUT,"| gzip > $out.gz" || die $!;
}

$len ||= 150;

#initialization of @Median
my @Median = ();
for my $i(0..($len-1))
{
	$Median[$i] = [];
}

sub Find_Median
{
	my $remainder = $#_%2;
	my $quotient = int($#_/2);
	my $the_Median;
	if($remainder == 0)
	{
		$the_Median = $_[$quotient];
	}
	else
	{
		$the_Median = int( ( $_[$quotient+1] + $_[$quotient] ) / 2 );
	}
	return $the_Median;
}

sub Out_Queue
{
	my $median_temp = shift @_;
	my $n = scalar(@$median_temp);
	my @out;
	for(my $j = 0; $j <= $n - 1 ; $j ++)
	{	
		#print "$median_temp->[$j]\t"; ## @$median_temp[$j]
		last if(!$median_temp->[$j]);
		push(@out, $median_temp->[$j]);
	}
	#print "$n\n";
	return @out;
}

my $gap = 0;
my $aloci = 0;
my @median = ();
my $median_value = 0;
while(<IN>)
{
	chomp;
	my @temp = split /\t/, $_;
	next if($temp[2] ne $chr);
	$gap = $temp[3] - $aloci;
	if($aloci != 0 and $gap != 0)
	{
		#progress of median finding
		for(my $i = 1; $i <= $gap; $i ++)
		{
			if(defined $Median[0]->[0]){
				@median = &Out_Queue(@Median);
				shift @Median;
				@median = sort{ $a<=>$b } @median;
				$median_value = &Find_Median(@median);
				my $print_loci = $aloci + $i - 1;
				print OUT "$print_loci\t$median_value\n";
				@median = ();
				$Median[$len-1] = [];
			}
			else{
				shift @Median;
				@median = ();
				$Median[$len-1] = [];
			}

			last if(!(defined $Median[0]->[0]));
		}
	}
	my @quality = split '', $temp[10];
	
	#adding the new qualities
	for my $k(0..$#quality)
	{
		my $quality_score = ord($quality[$k]) - 33;
		push(@{$Median[$k]}, $quality_score);
	}
	$aloci = $temp[3];
}

#print out the rest elements in @Median
while($Median[0]->[0]){
	@median = &Out_Queue(@Median);
	shift @Median;
	@median = sort {$a<=>$b} @median;
	$median_value = &Find_Median(@median);
	print OUT "$aloci\t$median_value\n";
	$aloci++;
}

close IN;
close OUT;
