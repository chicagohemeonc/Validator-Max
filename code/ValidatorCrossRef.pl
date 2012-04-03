#!/usr/bin/perl

# Written by: Jonathan Goya
# This script filters CPAS Xpress peptide reports with Sam Volchenboum's Validator pair-wise peptide identification software output, reporting only peptides of validated identity.

use warnings;
use strict;
($ARGV[0] ne "") || die "ValidatorCrossRef.pl validator_output.txt CPAS_peptides.txt crossref_pass.txt [Val1 | Val3]\n";


open(TRANSLATE1, "| tr '\r' '\n' < $ARGV[0] > $ARGV[0]_unix");
close TRANSLATE1;

open(TRANSLATE2, "| tr '\r' '\n' < $ARGV[1] > $ARGV[1]_unix");
close TRANSLATE2;

open(VALIDATOR,"$ARGV[0]_unix");
my %validatorhash;
my $valscan1index;
my $valscan2index;
my $valseq1index;
my $valseq2index;
my $pairindex;
my $mcatindex;

my $valheaderline = <VALIDATOR>;

chomp($valheaderline);
#print $valheaderline;
print "Parsing Validator output ...";	
my @valheader = split(/,/,$valheaderline);

if ($ARGV[3] eq "Val3")
        {	
        #print "########\n";
        #print join ("\n",@valheader);
        for (my $index = 0; $index <=$#valheader; $index ++)
	        {
	        $valheader[$index] =~ s/\s*//g;
	        if ($valheader[$index] eq "Pair")
		        {$pairindex = $index}
	        if ($valheader[$index] eq "ScanL")
		        {$valscan1index = $index}
	        if ($valheader[$index] eq "ScanH")
		        {$valscan2index = $index}
	        if ($valheader[$index] eq "PeptideL")
		        {$valseq1index = $index}
	        if ($valheader[$index] eq "PeptideH")
		        {$valseq2index = $index}
	        if ($valheader[$index] eq "MatchCategory")
		        {$mcatindex = $index}
	        }
#        print $valscan1index."\t".$valscan2index."\t".$valseq1index."\n";
        while (my $validatorline = <VALIDATOR>)
	        {
	        chomp($validatorline);
	        my @validator = split(/\t/,$validatorline);
	        my $pair = $validator[$pairindex];
	        my $scan1 = $validator[$valscan1index];
	        my $scan2 = $validator[$valscan2index];
	        my $valseq1 = $validator[$valseq1index];
	        my $valseq2 = $validator[$valseq2index];
	        my $mcat = $validator[$mcatindex];
	        if ($mcat eq "T")
	                {
        #	        print $valseq."\t".$scan1."\t".$scan2."\n";
	                $validatorhash{$pair}->{$scan1}->{'seq'}=$valseq1;
	                $validatorhash{$pair}->{$scan2}->{'seq'}=$valseq2;
	                }
	        elsif ($mcat eq "L")
	                {
        #	        print $valseq."\t".$scan1."\t".$scan2."\n";
	                $validatorhash{$pair}->{$scan1}->{'seq'}=$valseq1;
	                }
                elsif ($mcat eq "H")
	                {
        #	        print $valseq."\t".$scan1."\t".$scan2."\n";
	                $validatorhash{$pair}->{$scan2}->{'seq'}=$valseq2;
	                }
	        }
        }
elsif ($ARGV[3] eq "Val1")
        {
	
        #print "########\n";
        #print join ("\n",@valheader);

        for (my $index = 0; $index <=$#valheader; $index ++)
	        {
	        $valheader[$index] =~ s/\s*//g;
			#print $index,$valheader[$index]."\n";
	        if ($valheader[$index] eq "Scan1")
		        {$valscan1index = $index}
	        if ($valheader[$index] eq "Scan2")
		        {$valscan2index = $index}
	        if ($valheader[$index] eq "Peptide")
		        {$valseq1index = $index}
	        }
        #print $valscan1index."\t".$valscan2index."\t".$valseq1index."\n";
		#exit;
        my $pair =1;
        while (my $validatorline = <VALIDATOR>)
	        {
	        chomp($validatorline);
	        my @validator = split(/,/,$validatorline);
	        my $scan1 = $validator[$valscan1index];
	        my $scan2 = $validator[$valscan2index];
	        my $valseq1 = $validator[$valseq1index];
                $validatorhash{$pair}->{$scan1}->{'seq'}=$valseq1;
                $validatorhash{$pair}->{$scan2}->{'seq'}=$valseq1;
                $pair ++;
		
	        }
        
        
        
        }
else
        {
        die "Please indicate whether you are attempting to parse a Validator 1 or Validator 3 file: \nValidatorCrossRef.pl validator_output.txt CPAS_peptides.txt crossref_pass.txt Val[1|3]\n";
        }
print "done.\n";
print "Cross-referencing CPAS output ...";


open(PEPREP,"$ARGV[1]_unix");					#open Scaffold Spectrum Report
open(PEPOUT,">$ARGV[2]");
#open(PEPFAIL,">$ARGV[3]");
my $scanindex;
my $chargeindex;
my $seqindex;

while (my $headerline = <PEPREP>)
        {
        print PEPOUT $headerline;
        #print PEPFAIL $headerline;
        chomp($headerline);
        $headerline =~ tr/\"//d;
        my @header = split(/\t/,$headerline);	
        if ($header[0] and $header[0] eq "Scan")		#find header line
	        {
	        #print "########\n";
	        for (my $index = 0; $index <=$#header; $index ++)
		        {
		        if ($header[$index] eq "Scan")
			        {$scanindex = $index}
		        if ($header[$index] eq "Z")
			        {$chargeindex = $index}
		        if ($header[$index] eq "Peptide")
			        {$seqindex = $index}
		        }
			
	        while (my $pepline = <PEPREP>)
		        {
		        chomp($pepline);		#parse scan ID lines
		        my $peplineraw = $pepline;		
		        $pepline =~ tr/\"//d;
		        my @pep = split(/\t/,$pepline);
	
		        if ($pep[0] >= 1)
			        {
			        my $charge = $pep[$chargeindex];
			        my $scan = $pep[$scanindex];
			        my $seqfull = $pep[$seqindex];
			        my @seq = split(/\./,$seqfull);
			        my $sequence = $seq[1];
			        $sequence =~ tr/[\^|\"|\#|\~|\']//d;
			        my @pairs = sort keys %validatorhash;
			        foreach my $pair (@pairs)
			                {
		                        if ($validatorhash{$pair}->{$scan} and $validatorhash{$pair}->{$scan}->{'seq'} eq $sequence)
				                {
				                $validatorhash{$pair}->{$scan}->{'CPASline'}=$peplineraw;
				                }
				        }
			        }
		        }
		}
	}
my %crossref;
my @pairs = sort keys %validatorhash;
foreach my $pair (@pairs)
        {
        my @scans = sort keys %{$validatorhash{$pair}};
        if (scalar(@scans) == 2 && 
            $validatorhash{$pair}->{$scans[0]}->{'CPASline'} &&
            $validatorhash{$pair}->{$scans[1]}->{'CPASline'})
	        {
	        $crossref{$scans[0]}->{'CPASline'}=$validatorhash{$pair}->{$scans[0]}->{'CPASline'};
	        $crossref{$scans[1]}->{'CPASline'}=$validatorhash{$pair}->{$scans[1]}->{'CPASline'};
	        }
	if (scalar(@scans) == 1 && 
            $validatorhash{$pair}->{$scans[0]}->{'CPASline'})
	        {
	        $crossref{$scans[0]}->{'CPASline'}=$validatorhash{$pair}->{$scans[0]}->{'CPASline'};
	        }	
        }

        
my @scans = sort { $a <=> $b }keys %crossref;
foreach my $scan (@scans)
        {
        print PEPOUT $crossref{$scan}->{'CPASline'}."\n";
        
        }
print "done.\n";
unlink ("$ARGV[0]_unix");
unlink ("$ARGV[1]_unix");

				
