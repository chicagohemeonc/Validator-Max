#!/usr/bin/perl -w

use strict;
use Spreadsheet::ParseExcel;
my @filelist = "";
my $filename = "";
use FindBin qw($Bin); 
open (FILE, $ARGV[0]);
@filelist = <FILE>;
close(FILE);
my $data_path = $filelist[0];
chomp($data_path);
my $path = $Bin . "/../output_files/" . $data_path;


foreach $filename (@filelist[1..$#filelist])
{

chomp $filename;
print "Analyzing: $filename ...\n";
my ($iPepH, $iPepL, $iZM, $iIsoModsM, $iMatch, $iBYScore, $iBYScoreLL, $iBYScoreHH, $iProtL, $iProtH);

my $oExcel = new Spreadsheet::ParseExcel;

$filename =~ /(.+).xls/;
my $file_core = $1;
my $filename2="$path"."$file_core"."_MascotMatch.tsv";
my $filename3="$path"."$file_core"."_HLRecover.tsv";
my $filename4="$path"."$file_core"."_AllMatch.tsv";
open (MASCOTMATCH, ">$filename2");
open (HLRECOVER, ">$filename3");
open (ALLMATCH, ">$filename4");

my $oBook = $oExcel->Parse($path.$filename);


my($iR, $iC, $oWkS, $oWkC);
print "FILE  :", $oBook->{File} , "\n";
print "COUNT :", $oBook->{SheetCount} , "\n";
#print "AUTHOR:", $oBook->{Author} , "\n"
 #if defined $oBook->{Author};

for(my $iSheet=0; $iSheet < $oBook->{SheetCount} ; $iSheet++)
{
 $oWkS = $oBook->{Worksheet}[$iSheet];
 #print $oWkS->{Name}."\n";
 if ($oWkS->{Name} eq "pair_data")
 {
 #print "Sheet Name TRUE\n";
  for(my $iCh = $oWkS->{MinCol} ;
      defined $oWkS->{MaxCol} && $iCh <= $oWkS->{MaxCol} ;
      $iCh++)  
  {
   my $header = $oWkS->{Cells}[0][$iCh]->unformatted;
   print MASCOTMATCH $header."\t";
   print HLRECOVER $header."\t";
   print ALLMATCH $header."\t";
   #print $iCh."\t".$header."\n";
   if ($oWkS->{Cells}[0][$iCh]->unformatted eq "Charges match")
   {$iZM=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq "Peptide L")
   {$iPepL=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq "Peptide H")
   {$iPepH=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq 'Match?')
   {$iMatch=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq 'Isotope + Mods match L')
   {$iIsoModsM=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq "Mascot proteins L")
   {$iProtL=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq "Mascot proteins H")
   {$iProtH=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq 'B/Y Match Score')
   {$iBYScore=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq 'B/Y Match Score L L')
   {$iBYScoreLL=$iCh}
   elsif ($oWkS->{Cells}[0][$iCh]->unformatted eq 'B/Y Match Score H H')
   {$iBYScoreHH=$iCh}
  }
  #print join ("\t", ($iPepH, $iPepL, $iZM, $iIsoModsM, $iMatch, $iBYScore, $iBYScoreLL, $iBYScoreHH))."\n";
  print MASCOTMATCH "Match Category\n";
  print HLRECOVER "Match Category\n";
  print ALLMATCH "Match Category\n";
  for(my $iR = 1 ;
      defined $oWkS->{MaxRow} && $iR <= $oWkS->{MaxRow} ;
      $iR++)
  {
   if (($oWkS->{Cells}[$iR][$iMatch]->unformatted eq 'T') &&
       ($oWkS->{Cells}[$iR][$iZM]->unformatted eq 'T') &&
       ($oWkS->{Cells}[$iR][$iIsoModsM]->unformatted eq 'T') &&
       ($oWkS->{Cells}[$iR][$iBYScore]->unformatted >= 10))
   {
    for(my $iC = $oWkS->{MinCol} ;
        defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
        $iC++)
    {
     $oWkC = $oWkS->{Cells}[$iR][$iC];
     print MASCOTMATCH $oWkC->unformatted."\t";
     print ALLMATCH $oWkC->unformatted."\t";
    }
    print MASCOTMATCH "T\n";
    print ALLMATCH "T\n";
   }
   
   if ($oWkS->{Cells}[$iR][$iZM]->unformatted eq 'T')
   {
   unless ($oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted and
           $oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted and 
           $oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted >= 15 and
           $oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted >= 15)
   {
    if ($oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted and
        $oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted >= 15)
    {
     for(my $iC = $oWkS->{MinCol} ;
         defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
         $iC++)
     { 
      $oWkC = $oWkS->{Cells}[$iR][$iC];
      if ($iC == $iPepH)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iPepL]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iPepL]->unformatted."\t";
      }
      elsif ($iC == $iProtH)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iProtL]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iProtL]->unformatted."\t";
      }
      else
      {
       print HLRECOVER $oWkC->unformatted."\t";
       print ALLMATCH $oWkC->unformatted."\t";
      }
     }
     print HLRECOVER "L\n";
     print ALLMATCH "L\n";
    }

    if ($oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted and 
        $oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted >= 15)
    {
     for(my $iC = $oWkS->{MinCol} ;
         defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
         $iC++)
     {
      $oWkC = $oWkS->{Cells}[$iR][$iC];
      if ($iC == $iPepL)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iPepH]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iPepH]->unformatted."\t";
      }
      elsif ($iC == $iProtL)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iProtH]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iProtH]->unformatted."\t";
      }
      else
      {
       print HLRECOVER $oWkC->unformatted."\t";
       print ALLMATCH $oWkC->unformatted."\t";
      }
     }
     print HLRECOVER "H\n";
     print ALLMATCH "H\n";
    }
   } 
   elsif ($oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted > $oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted)
   {
    {
     for(my $iC = $oWkS->{MinCol} ;
         defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
         $iC++)
     { 
      $oWkC = $oWkS->{Cells}[$iR][$iC];
      if ($iC == $iPepH)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iPepL]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iPepL]->unformatted."\t";
      }
      elsif ($iC == $iProtH)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iProtL]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iProtL]->unformatted."\t";
      }
      else
      {
       print HLRECOVER $oWkC->unformatted."\t";
       print ALLMATCH $oWkC->unformatted."\t";
      }
     }
     print HLRECOVER "L\n";
     print ALLMATCH "L\n";
    }
   }
   elsif ($oWkS->{Cells}[$iR][$iBYScoreHH]->unformatted > $oWkS->{Cells}[$iR][$iBYScoreLL]->unformatted)
   {
    {
     for(my $iC = $oWkS->{MinCol} ;
         defined $oWkS->{MaxCol} && $iC <= $oWkS->{MaxCol} ;
         $iC++)
     {
      $oWkC = $oWkS->{Cells}[$iR][$iC];
      if ($iC == $iPepL)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iPepH]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iPepH]->unformatted."\t";
      }
      elsif ($iC == $iProtL)
      {
       print HLRECOVER $oWkS->{Cells}[$iR][$iProtH]->unformatted."\t";
       print ALLMATCH $oWkS->{Cells}[$iR][$iProtH]->unformatted."\t";
      }
      else
      {
       print HLRECOVER $oWkC->unformatted."\t";
       print ALLMATCH $oWkC->unformatted."\t";
      }
     }
     print HLRECOVER "H\n";
     print ALLMATCH "H\n";
    }
   }
  }
  }
 }
}
print "done\n";
}
