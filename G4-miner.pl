#!/usr/bin/perl -w
#########################################################################
# FileName: G4-miner.pl
# Version: 49916fbf-51f7-4d04-abd1-42b239642890
# Author: Tu's Lab <jtu@seu.edu.cn>
# CreatedTime: Tue 04 Sep 2021 04:20:44 PM CST
#########################################################################
use strict;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $red = "\033[0;31m";
my $end = "\033[0m";

my $homepath = dirname(File::Spec->rel2abs( $0 ));
$homepath=~s/^\s+|\s+$//g;

my $scriptpath = "$homepath/Scripts";
my $path1 = "$scriptpath/1.preprocess";
my $path2 = "$scriptpath/2.foward";
my $path3 = "$scriptpath/3.reverse";

my ($in, $ref, @chr, $clean, $hdir, $read_length, $region_length, $specie, $samp, $py2, $G4Predict, $help);

GetOptions
(
	"i|in=s"=>\$in,
	"d|dir=s"=>\$hdir,
	"g|ref=s"=>\$ref,
	"n|samp=s"=>\$samp,
	"rm=i"=>\$clean,
	"ch=s"=>\@chr,
	"py=s"=>\$py2,
	"g4=s"=>\$G4Predict,
	"rl=i"=>\$read_length,
	"dl=i"=>\$region_length,
	"sp=s"=>\$specie,
	"help|?"=>\$help,
);

my $usage=<<INFO;
Usage:
	perl $0 [options]
Options:
		-i <file> <${red}input bam/sam file$end>
		-d <path> <${red}the root directory for saving output$end>
		-g <file> <${red}the fasta file of the whole genome reference$end>
		-rm <int> <${red}rm the directory of \$hdir/\$samp, default: 1, yes$end>
		-ch <string> <${red}specify some specific chromosomes for analysing, -ch chr1 -ch chr2 ... -ch chrX> <Optional, Default genome wide$end>
		-py <file> <the absolute path of python2> (default: `which python2`)
		-g4 <file> <the absolute path of g4predict> (default: `which g4predict`)
		-rl <int> <the longest length of reads, ${red}default = 150 nt$end>
		-dl <int> <the length of detection region, ${red}defalut = 75 nt$end>
		-sp <float> <${red}Specie for selecting cutoff of G ratio, here have six species: human/mouse/fruitfly/thaliana/elegans/yeast,
		             default: human with 0.28, if the specie did not include in this paper, can calculate by one more step$end>
		-n <string> <${red}the name of the sample$end>
INFO

die $usage if ($help || !$in || !$hdir || !$samp || !$ref);

$read_length ||= 150;

$region_length ||= 75;

if(!(defined($clean)))
{
	$clean = 1;
}

$specie ||= "human";

die "${red}$in !exists\n" if(!(-e "$in"));

die "${red}$ref !exists\n" if(!(-e "$ref"));

my %gratio = (-human=>0.28, -mouse=>0.28, -fruitfly=>0.27, -thaliana=>0.24, -elegans=>0.25, -yeast=>0.29,);

my %cratio = (-human=>0.28, -mouse=>0.28, -fruitfly=>0.27, -thaliana=>0.24, -elegans=>0.25, -yeast=>0.29,);

my $cmd;

print $red,"## Split reference fa to chromosome fa && run g4predict$end\n";
$py2 ||= `which python2`;
$py2=~s/^\s+|\s+$//g;
die "${red}There did not install python2 $py2$end\n" if(!$py2);

$G4Predict ||= `which g4predict`;
$G4Predict=~s/^\s+|\s+$//g;
die "${red}There did not install g4predict $G4Predict$end\n" if(!$G4Predict);

if($ref =~ /\.gz$/)
{
	open IF, "gzip -dc $ref |" || die $!;
}
else
{
	open IF,"$ref" || die $!;
}

my @chrs = ();

my $odir = "$hdir/$samp";

`rm -rf $hdir/$samp/*` if(-e "$odir" && $clean);

if(! (-e "$odir") )
{
	`mkdir -p $odir` and die "${red}Error: can't create $odir$end\n";
}

if(@chr)
{
	@chrs = @chr;
}

my $tag = 0;
while (<IF>)
{
	chomp;
	my $chr;
	my $out;
	if ($_ =~ />(.*)/)
	{
		$chr = $1;
		push(@chrs, $chr) if(!@chr);
		$tag = 0;
		$out = "$chr.fa";
		if(!(-e "$odir/$chr/$out"))
		{
			if(!@chr)
			{
				$tag ++;
				if(! (-e "$odir/$chr") )
				{
					`mkdir -p $odir/$chr` and die "${red}Error: can't create $odir/$chr$end\n";
				}
				if($clean)
				{
					open OUT, "> $odir/$chr/$out" || die $!;
					print OUT "$_\n";
				}
			}
			else
			{
				my $grep = grep { $_ eq $chr } @chr;
				$grep=~s/^\s+|\s+$//g;
				if($grep){
					$tag ++;
					if(! (-e "$odir/$chr") )
					{
						`mkdir -p $odir/$chr` and die "${red}Error: can't create $odir/$chr$end\n";
					}
					if($clean)
					{
						open OUT, "> $odir/$chr/$out" || die $!;
						print OUT "$_\n";
					}
				}
				else
				{
					$tag = 0;
				}
			}
		}
		else
		{
			$tag = 0;
		}
	}
	else
	{
		if($clean)
		{
			print OUT "$_\n" if $tag;
		}
	}
}
close IF;
close OUT;

print $red,"## Judge the fasta index existx or not$end\n";
my $fai = "$ref.fai";

my $fai2 = "${ref}i";

if(!(-e "$fai"))
{
	if(!(-e "$fai2"))
	{
		$cmd = "samtools faidx $ref";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
		$fai = "$ref.fai";
	}
	else
	{
		$fai = $fai2;
	}
}

if($clean)
{
	print $red,"## Generate reverse and forward reads and save into dirrerent file$end\n";
	$cmd = "perl $scriptpath/ForRev.pl -i $in -f $fai -or $odir/$samp.r.sort.bam -of $odir/$samp.f.sort.bam";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";


	for my $chr(@chrs)
	{
		print $red,"## Predict G4 for each chromosome by using $G4Predict$end\n";
		my $fa = "$odir/$chr/$chr.fa";
		$cmd = "$py2 $G4Predict intra -f $fa -b $odir/$chr/PQ_g4predict.bed -s -M";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
	}
}

my $g;

my $c;

if(exists $gratio{-$specie})
{
	$g = $gratio{-$specie};
	
	$c = $cratio{-$specie};
}
else
{
	$g = &PQ_Gratio($odir, @chrs);
	
	$c = &PQ_Cratio($odir, \@chrs);
}

print "PQ_Gratio is $g\n";

print "PQ_Cratio is $c\n";

sub PQ_Gratio
{
	my @lists = @_;
	my $dir_main = $lists[0];
	my @num;
	my $sum_pq = 0;
	my $tag = -1;

	for my $chr (@lists)
	{
		$tag ++;
		next unless($tag >= 1);
		open IN,"$dir_main/$chr/PQ_g4predict_plus" || die $!;
		open INF,"$dir_main/$chr/${chr}.fa" || die $!;

		my @ref = ('N');
		while(<INF>)
		{
			next if $_ =~ />/;
			chomp;
			my @seq = split("",$_);
			push(@ref, @seq);
		}
		close INF;

		while(<IN>)
		{
			chomp;
			$sum_pq ++;
			my @temp = split/\s+/, $_;
			my $pq_s = $temp[1] - 62;
			my $pq_e = $temp[1] + 12;
			my $count = 0;
			for(my $i = $pq_s; $i <= $pq_e; $i++)
			{
				if(!$ref[$i])
				{
					last;
				}
				
				if($ref[$i] eq "G")
				{
					$count ++;
					$num[$count] ++;
				}
			}
		}
		close IN;
	}

	my $per99 = $sum_pq * 0.99;
	my $gnum = 0;
	for my $i (0..$#num)
	{
		last if($num[$i] < $per99);
		$gnum ++ if $num[$i] >= $per99; 
	}
	
	my $gratio = sprintf "%.2f", $gnum / 75;

	return $gratio;
}

sub PQ_Cratio
{
    my ($dir_main, $chrs) = @_; 
    my @num;
    my $sum_pq = 0;

    for my $chr (@$chrs)
    {
        open IN,"$dir_main/$chr/PQ_g4predict_minus" || die $!;
        open INF,"$dir_main/$chr/${chr}.fa" || die $!;

        my @ref = ('N');
        while(<INF>)
        {
            next if $_ =~ />/;
            chomp;
            my @seq = split("",$_);
            push(@ref, @seq);
        }
        close INF;

        while(<IN>)
        {
            chomp;
            $sum_pq ++;
            my @temp = split/\s+/, $_;
            my $pq_s = $temp[0] - 12;
            my $pq_e = $temp[0] + 64;
            my $count = 0;
            for(my $i = $pq_s; $i <= $pq_e; $i++)
            {
                if(!$ref[$i])
                {
                    last;
                }

                if($ref[$i] eq "C")
                {
                    $count ++;
                    $num[$count] ++;
                }
            }
        }
        close IN;
    }

    my $per99 = $sum_pq * 0.99;
    my $cnum = 0;
    for my $i (0..$#num)
    {
        last if($num[$i] < $per99);
        $cnum ++ if $num[$i] >= $per99;
    }

    my $cratio = sprintf "%.2f", $cnum / 75;

    return $cratio;
}

my $tag_step = 0;

if($clean)
{
	for my $chr(@chrs)
	{	
		$tag_step ++;
		my $fa = "$odir/$chr/$chr.fa";
		$cmd = "perl $path1/PQ_g4predict_seq.pl -i $odir/$chr/PQ_g4predict.bed -f $fa -o PQ_g4predict";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		print $red,"## Count the median of the quality value of each base on $chr with $odir/$samp.r.sort.bam$end\n" if($tag_step == 1);
		print $red,"## Start Analysis read mapped to reverse strand of genome$end\n" if($tag_step == 1);
		$cmd = "perl $path1/Median.pl -i $odir/$samp.r.sort.bam -o $odir/$chr/Medianr.gz -l $read_length -ch $chr";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path2/Positive_region_median.pl -i $odir/$chr/PQ_g4predict_plus -m $odir/$chr/Medianr.gz -o $odir/$chr/positive_array_r -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path2/Negative_sloci.pl -i $odir/$chr/PQ_g4predict_plus -m $odir/$chr/Medianr.gz -o $odir/$chr/negative_sloci_r.gz -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path2/Negative_seq_filter.pl -i $odir/$chr/negative_sloci_r.gz -f $fa -o $odir/$chr/false_positive_sloci_r.gz -c $g";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path2/Negative_count.pl -i $odir/$chr/negative_sloci_r.gz -m $odir/$chr/Medianr.gz -o $odir/$chr/negative_array_r -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path2/Negative_count.pl -i $odir/$chr/false_positive_sloci_r.gz -m $odir/$chr/Medianr.gz -o $odir/$chr/false_positive_array_r -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		print $red,"## Count the median of the quality value of each base on $chr with $odir/$samp.f.sort.bam$end\n" if($tag_step == 1);
		print $red,"## Start Analysis read mapped to forward strand of genome$end\n" if($tag_step == 1);
		$cmd = "perl $path1/Median.pl -i $odir/$samp.f.sort.bam -o $odir/$chr/Medianf.gz -l $read_length -ch $chr";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path3/Positive_region_median.pl -i $odir/$chr/PQ_g4predict_minus -m $odir/$chr/Medianf.gz -o $odir/$chr/positive_array_f -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path3/Negative_sloci.pl -i $odir/$chr/PQ_g4predict_minus -m $odir/$chr/Medianf.gz -o $odir/$chr/negative_sloci_f.gz -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path3/Negative_seq_filter.pl -i $odir/$chr/negative_sloci_f.gz -f $fa -o $odir/$chr/false_positive_sloci_f.gz -c $c -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path3/Negative_count.pl -i $odir/$chr/negative_sloci_f.gz -m $odir/$chr/Medianf.gz -o $odir/$chr/negative_array_f -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

		$cmd = "perl $path3/Negative_count.pl -i $odir/$chr/false_positive_sloci_f.gz -m $odir/$chr/Medianf.gz -o $odir/$chr/false_positive_array_f -l $region_length";
		print "\nNOTICE: Running with system command <$cmd>\n";
		system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
	}
}

if($clean)
{
	for my $pre("false_positive_array", "positive_array", "negative_array")
	{
		if(!@chrs)
		{
			$cmd = "perl $path2/array_merged.pl -i ${pre}_r -g $ref -d $odir";
			print "\nNOTICE: Running with system command <$cmd>\n";
			system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

			$cmd = "perl $path3/array_merged.pl -i ${pre}_f -g $ref -d $odir";
			print "\nNOTICE: Running with system command <$cmd>\n";
			system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
		}
		else
		{
			my $chrs_used;
			for my $i(@chrs)
			{
				$chrs_used .= (!$chrs_used) ? "$i" : ",$i";
			}

			$cmd = "perl $path2/array_merged.pl -i ${pre}_r -cs $chrs_used -d $odir";
			print "\nNOTICE: Running with system command <$cmd>\n";
			system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

			$cmd = "perl $path3/array_merged.pl -i ${pre}_f -cs $chrs_used -d $odir";
			print "\nNOTICE: Running with system command <$cmd>\n";
			system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
		}
	}
}

my $gzr;
my $gzf;
my $pqr;
my $pqf;

for my $chr(@chrs)
{
	my $dirss = "$odir/$chr";
	$gzr .= (!$gzr) ? "$dirss/negative_sloci_r.gz" : " $dirss/negative_sloci_r.gz";
	$gzf .= (!$gzf) ? "$dirss/negative_sloci_f.gz" : " $dirss/negative_sloci_f.gz";
	$pqr .= (!$pqr) ? "$dirss/PQ_g4predict_plus" : " $dirss/PQ_g4predict_plus";
	$pqf .= (!$pqf) ? "$dirss/PQ_g4predict_minus" : " $dirss/PQ_g4predict_minus";
}

my $lcgzr = `gzip -dc $gzr | wc -l`;
$lcgzr=~s/^\s+|\s+$//g;

my $lcgzf = `gzip -dc $gzf | wc -l`;
$lcgzf=~s/^\s+|\s+$//g;

my $lcpqr = `cat $pqr | wc -l`;
$lcpqr=~s/^\s+|\s+$//g;

my $lcpqf = `cat $pqf | wc -l`;
$lcpqf=~s/^\s+|\s+$//g;

my @fplr = &CoutArry($lcgzr, "$odir/false_positive_array_r");
my @fplf = &CoutArry($lcgzf, "$odir/false_positive_array_f");
my @plr = &CoutArry($lcpqr, "$odir/positive_array_r");
my @plf = &CoutArry($lcpqf, "$odir/positive_array_f");

my ($mr, $nr, $maxr) = &MaxMN(\@fplr, \@plr);
my ($mf, $nf, $maxf) = &MaxMN(\@fplf, \@plf);

sub MaxMN
{
	my $m = 12;
	my $n = 20;
	my ($fpl, $pl) = @_;

	for(my $i = 1; $i <= $n; $i++)
	{
		for(my $j = 1; $j <= $m; $j++)
		{
			$pl->[$i]->[$j] = 0 if($fpl->[$i]->[$j] >= 0.01);
		}
	}

	my ($mx, $nx);
	my $max = 0;
	for(my $i = 1; $i <= $n; $i++)
	{
		for(my $j = 1; $j <= $m; $j++)
		{
			if($pl->[$i]->[$j] > $max)
			{
				$mx = $j;
				$nx = $i;
				$max = $pl->[$i]->[$j];
			}
		}
	}

	return ($mx, $nx, $max);
}

sub CoutArry
{
	my ($lc, $file) = @_;
	my @array;

	open IN,"$file" || die "cannot open this file!\n";
	my $row = 1;
	while(<IN>)
	{
		chomp;
		my @temp = split/\s+/, $_;
		for my $j (0..$#temp)
		{
			$array[$row]->[$j+1] = $temp[$j]/$lc;
		}
		$row ++;
	}
	close IN;

	return @array;
}

for my $chr(@chrs)
{
	my $fa = "$odir/$chr/$chr.fa";
	my $dir = "$odir/$chr";
	$cmd = "perl $path2/PeaksExtract.pl -i $dir/Medianr.gz -f $fa -c $mr -q $nr -o $dir/OQ_r.gz";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path2/G_ratio_cutoff.pl -i $dir/OQ_r.gz -c $g -o $dir/G4_OQ_unmerged_r";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path2/merge.pl -i $dir/G4_OQ_unmerged_r -f $fa -o $dir/G4_OQ_merged_r";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path2/PQ_in_OQ_g4predict.pl -1 $dir/PQ_g4predict_plus -2 $dir/G4_OQ_merged_r -o $dir/G4_PQinOQ_r";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path3/PeaksExtract.pl -i $dir/Medianf.gz -f $fa -c $mf -q $nf -o $dir/OQ_f.gz";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path3/C_ratio_cutoff.pl -i $dir/OQ_f.gz -c $c -o $dir/G4_OQ_unmerged_f";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path3/merge.pl -i $dir/G4_OQ_unmerged_f -f $fa -o $dir/G4_OQ_merged_f";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";

	$cmd = "perl $path3/PQ_in_OQ_g4predict.pl -1 $dir/PQ_g4predict_minus -2 $dir/G4_OQ_merged_f -o $dir/G4_PQinOQ_f";
	print "\nNOTICE: Running with system command <$cmd>\n";
	system ($cmd) and die $red,"Error running system command: <$cmd>$end\n";
}

print $red."Finished $samp G4-miner$end\n";

exit 0;