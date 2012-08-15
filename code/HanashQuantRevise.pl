####
#Copyright 2011 Samuel Volchenboum, Jonathan Goya, Gene Selkov, Chaim Kirby, 
#
#This file is part of Validator.
#
#Validator is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Validator is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with Validator.  If not, see <http://www.gnu.org/licenses/>.
#####
#!/usr/bin/perl -w
use strict;
use Cwd;

my %Normaliz; 
my $split;
# statistics sub_routine
use FindBin qw($Bin); 

my $dir = $Bin . "/supporting_data";


sub median {
    my @pole = @_;

    my $ret;

    @pole = sort(@pole);

    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }

    return $ret;
}



sub stat
	{
	my ($key, @ratios) = @_;
#	print $key."\n";
	my $normaliz;
	if ($split eq "all")
	        {$normaliz = $Normaliz{'Median'};}
        elsif ($split eq "sec")
                {
                my ($protein, $section) = split (/\t/, $key);
	        $normaliz = $Normaliz{$section}->{'Median'};
	        }
	
	my $total=0;
	my $average=0;
	my $stdev=0;
	my $sum_std=0;
	my $t_test=0;
	my $P_dF;
	my $P_value=1;
									#LINE_30
	my $c=0;
	my $null_hyp=0;
	
	if (scalar(@ratios) > 1) 
		{
		foreach my $val (@ratios) 
			{
			$val=$val/$normaliz;
			#print "val "."$val"."\n";
			$total=$total+log($val);
			#print "total "."$total"."\n";
								#LINE_40
			}
		$average=($total/scalar(@ratios));
		#print "Average "."$average"."\n";

		foreach my $val (@ratios)
			{
			#print "@ratios";
			#print "val "."$val"."\n";
			$sum_std=$sum_std + ((log($val)-$average)*(log($val)-$average));
			#print "sum_std "."$sum_std"."\n";
								#LINE_50
			}

		$stdev=sqrt((1/(scalar(@ratios)-1)*$sum_std));
		#print "stdev "."$stdev"."\n";
		#print scalar(@ratios);
		if($stdev ne 0) {$t_test=($average-$null_hyp)/($stdev/sqrt(scalar(@ratios)))};
		#print "ttest "."$t_test"."\n";
		$average = exp($average);
		#print "average "."$average"."\n";

								#LINE_60
		open (T_TABLE, "<$dir/T_table.txt"); 
		while (my $line = <T_TABLE>)
			{
  			chomp $line;
			my @table = split ("\t", $line);
  			if(scalar(@ratios) > 120) {$P_value=0;}
  			else 
				{
				if ($t_test > $table[0] and $c==0){}
	  			elsif ($t_test < $table[0] and $c==0)
					{
	          			$P_value = $table[(scalar(@ratios)-1)];
	    				$c=1;
					}
	  			}
	  
	  		}
		close T_TABLE;
		}

	if (scalar(@ratios) ==1) 
		{
		$average = $ratios[0]/$normaliz;
		print 
		open (P_TABLE, "<$dir/P_table.txt"); 
		while (my $line = <P_TABLE>)
			{
	  		chomp $line;
	  		my @table = split ("\t", $line);
	  		if ($average > ($table[0]) and $c==0){}
	  		elsif ($average < $table[0] and $c==0)
				{
	    			$P_value = $table[1];
	    			$c=1;
				}
		  	}
	 	close P_TABLE;
		@ratios=(); 
	  	}



	#print "mean: ".exp($average)."\t";
	#print "STDEV: $stdev\t";
	#print "T-TEST: $t_test\n";
	#print "P-Value: $P_value\n";
	return ($average,$stdev,$t_test,$P_value);
	}

#end statistics_subroutine






my %Quant;
my %AllQuant;

$split = $ARGV[1];
#open(TRANSLATE1, "| tr -d '\r' < $ARGV[0] > unix_$ARGV[0]");
open(TRANSLATE1, "| tr '\r' '\n' < $ARGV[0] > $ARGV[0]_unix");
close TRANSLATE1;

open (FILE, "$ARGV[0]_unix");

my $protindex;
my $seqindex;
my $lightareaindex;
my $heavyareaindex;
my $modindex;
my $chargeindex;
my $mzxmlindex;


while (my $headerline = <FILE>) 
	{
	chomp($headerline);
	if ($headerline)
		{
		$headerline =~ tr/\"//d;
		print $headerline."\n";
					
		my @header = split(/\t/,$headerline);
	
		if (($header[0] eq "Protein") || ($header[0] eq "Gene Name") || ($header[0] eq "Scan"))		#find header line
			{
			for (my $index = 0; $index <=$#header; $index ++)

				{
				print $index."\t"."$header[$index]"."\n";
				if ($header[$index] eq "Fraction Name")
					{
					$mzxmlindex = $index;
					}
				if (($header[$index] eq "Gene Name"))
					{$protindex = $index}
				#elsif ($header[$index] eq "Protein")
        			#	{$protindex = $index}
				if (($header[$index] eq "Sequence") || ($header[$index] eq "Peptide"))
					{$seqindex = $index}
				if ($header[$index] eq "Charge" || ($header[$index] eq "Z"))
					{$chargeindex = $index}
				if ($header[$index] eq "UnlabeledArea")
					{$lightareaindex = $index}
				elsif ($header[$index] eq "Light Area")
					{$lightareaindex = $index}
				if ($header[$index] eq "Modifications")
					{$modindex = $index}
				if ($header[$index] eq "LabeledArea")
					{$heavyareaindex = $index}
				elsif ($header[$index] eq "Heavy Area")
					{$heavyareaindex = $index}
	#			if ($header[$index] eq '')
	#				{$index = $index}

				}
				
#			print "$split.$chargeindex\n";

#			print "$mzxmlindex.$protindex\n";
			while (my $line = <FILE>)
				{
				my $key;
				chomp($line);		#parse scan ID lines		
				$line =~ tr/\"//d;
				my @line = split(/\t/,$line);
			
				if ($line[$heavyareaindex] && $line[$lightareaindex] && ($line[$heavyareaindex] > 100) && ($line[$lightareaindex] > 100)) 
					{
#				        print "$split\n";
				        if ($split eq "all")
				                {
#				                print $line[$protindex];
				                $key = $line[$protindex];
				                }
				        elsif ($split eq "sec")
				                {
				                $key = $line[$protindex]."\t".$line[$mzxmlindex];
#				                print $key."\n";
				                }
					my $ratio = $line[$heavyareaindex]/$line[$lightareaindex];
					push (@{$Quant{$key}}, $ratio);
					if ($split eq "all")
				        {push (@{$Normaliz{'Quants'}}, $ratio);}
					elsif ($split eq "sec")
					{push (@{$Normaliz{$line[$mzxmlindex]}->{'Quants'}}, $ratio);}
				        }
				}
			}
		}
	}

$ARGV[0] =~ /(.+).(tsv|q6|txt|csv)/;
my $file_core = $1;
my $filename2=$file_core.".".$split.".events";
my $filename3=$file_core.".".$split.".quant";
open (OUTPUT, ">$filename2");
open (OUTFILE1, ">$filename3");
if ($split eq "all")
        {
        print OUTFILE1 "protein"."\t"."MEAN.HvL"."\t"."lnMEAN.HvL"."\t"."MEAN.LvH"."\t"."lnMEAN.LvH"."\t"."EVENTS"."\t"."STDEVln"."\t"."T-TESTln"."\t"."P-VALUEln"."\n";
        my @ratios = @{$Normaliz{'Quants'}};
        $Normaliz{'Median'} = median(@ratios);
        }
elsif ($split eq "sec")
        {
        print OUTFILE1 "protein"."\t"."FileSection"."\t"."MEAN.HvL"."\t"."lnMEAN.HvL"."\t"."MEAN.LvH"."\t"."lnMEAN.LvH"."\t"."EVENTS"."\t"."STDEVln"."\t"."T-TESTln"."\t"."P-VALUEln"."\n";
        foreach my $section (keys %Normaliz)
                {
                my @ratios = @{$Normaliz{$section}->{'Quants'}};
#                $Normaliz{$section}->{'Median'} = median(1,4,5,6,9,3,3);
                $Normaliz{$section}->{'Median'} = median(@ratios);
                print $section."\t".$Normaliz{$section}->{'Median'}."\n";
                }
        }


foreach my $output (keys %Quant)
	{
	print OUTPUT $output."\t@{$Quant{$output}}\n";
	my @results = &stat ($output, @{$Quant{$output}}); 				#send the data to the statistics subroutine
	if ($split eq "all")
	        {print OUTFILE1 $output."\t".$results[0]."\t".log($results[0])."\t".(1/$results[0])."\t".(-1*log($results[0]))."\t".scalar(@{$Quant{$output}})."\t".$results[1]."\t".$results[2]."\t".$results[3]."\n";} #result from stat (mean, stdev, t-test, pvalue)
        elsif ($split eq "sec")
		{print OUTFILE1 $output."\t".$results[0]."\t".log($results[0])."\t".(1/$results[0])."\t".(-1*log($results[0]))."\t".scalar(@{$Quant{$output}})."\t".$results[1]."\t".$results[2]."\t".$results[3]."\n";} #result from stat (mean, stdev, t-test, pvalue)
	}


close OUTFILE1;
close OUTPUT;
close FILE;
unlink "$ARGV[0]_unix";
