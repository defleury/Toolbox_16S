#!/usr/bin/perl
################################################################################
################################################################################
#Toolbox 16S: convenience script to process Infernal alignment outputs
#
#This script truncates an Infernal alignment in Pfam format to specified
#flanking positions and exports a Stockholm format alignment file, an alignment
#scores file and a list of removed sequences.
#
#Usage:
#filter.infernal.pl
#  <align.log>          #input: raw alignment log file
#  <align.tmp.sto>      #input: raw alignment sto file
#  <cut.start>          #parameter: alignment pruning start position (relative to Infernal model)
#  <cut.end>            #parameter: alignment pruning end position
#  <max.inserts>        #parameter: maximum number of inserts accepted
#  <max.N>              #parameter: maximum number of N characters accepted
#  <max.gaps>           #parameter: maximum number of gap characters accepted
#  <output.scores>      #output: alignment scores file
#  <output.removed>     #output: removed sequences
#  <output.sto>         #output: alignment in stockholm format
#
#sebastian.schmidt@embl.de 2016-12-07
################################################################################
################################################################################

use strict;

#Get input
my $file_align_log = shift @ARGV or die;
my $file_align_tmp = shift @ARGV or die;
my $cut_start = shift @ARGV or die;
my $cut_end = shift @ARGV or die;
my $max_inserts = shift @ARGV or die;
my $max_N = shift @ARGV or die;
my $max_gaps = shift @ARGV or die;
my $file_scores = shift @ARGV or die;
my $file_removed = shift @ARGV or die;
my $file_sto = shift @ARGV or die;

#Set parameters for flanking
my $start = $cut_start - 1;
my $length = $cut_end - $cut_start;

#Preallocate
my %Sequences; my %Filtered; my %Use;

#Process scores
open(RAW, $file_align_log) or die;
open(SCORES, ">$file_scores");
while (<RAW>) {
   chomp;
   #Edit whitspace (spaces to tabs)
   s/ +/\t/g;
   my @line = split /\t/;
   
   #Require positive alignment score against bacterial model
   if (($line[7] > 0) & ($line[12] > 0)) {
     print SCORES "$line[2]\tbacteria\t$line[7]\t$line[12]\n";
     $Sequences{$line[2]}{Score} = $line[7];
     $Use{$line[2]} = 1;
   }
   else {print SCORES "$line[2]\tnone\t0\t0\n"}
}
close RAW; close SCORES;

#Process sequences
open(SEQ, $file_align_tmp) or die;
open(ALIGNED, ">$file_sto") or die;
open(REMOVED, ">$file_removed") or die;
LINE:
while (<SEQ>) {
  chomp;
  #Edit whitspace (spaces to tabs)
  s/ +/\t/g;
  #Get current sequence
  my ($acc, $alignment) = split /\t/;
  my $tmp = $alignment;
  
  #Skip sequences with too many unaligned bases (i.e., >= max.inserts)
  if (($tmp =~ tr/acgun//) >= $max_inserts) {print REMOVED "$acc\tinsertions\t$tmp\n"; next LINE}
  
  #Remove lower case letters (insertions in sequence)
  $alignment =~ s/[acgun]//g;
  #Remove "points" (gaps in sequence)
  $alignment =~ s/\.//g;
  #Remove insertions
  $alignment =~ s/[a-z]//g;
  
  #Prune alignment
  my $pruned = substr $alignment, $start, $length;
  
  #Skip if sequence is too gappy
  $tmp = $pruned;
  if (($tmp =~ tr/-//) >= $max_gaps) {print REMOVED "$acc\tgaps\t$pruned\n"; next LINE}
  
  #Skip if sequence has two (or more) consecutive Ns
  if ($pruned =~ /N-*N/) {print REMOVED "$acc\tconsecutive_Ns\t$pruned\n"; next LINE}
  
  #Skip if sequence has too many Ns
  $tmp = $pruned;
  if (($tmp =~ tr/N//) >= $max_N) {print REMOVED "$acc\tN_count\t$pruned\n"; next LINE}
  
  print ALIGNED "$acc $pruned\n";
}
print ALIGNED "//\n";
