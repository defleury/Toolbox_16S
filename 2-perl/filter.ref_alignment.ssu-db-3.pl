#!/usr/bin/perl
################################################################################
################################################################################
#Toolbox 16S: convenience script to filter the full-ssu-db3 bacterial alignment
#down to the sequences used in mapref-db.
#
#Usage:
#./filter.re_alignment.ssu-db-3.pl <otu.file> <full.alignment> 
#
#sebastian.schmidt@embl.de 2016-12-21
################################################################################
################################################################################

use strict;
use warnings;
use Sort::Naturally;

#Get input
my $file_otu = shift @ARGV or die;
my $file_alignment_full = shift @ARGV or die;

#Parse full alignment
my %Seq;
open(FULL, $file_alignment_full) or die;
while (my $line = <FULL>) {
  chomp $line;
  if ($line =~ /^>/) {
    $line =~ s/^>//;
    my $seqline = <FULL>; chomp $seqline;
    $Seq{$line} = $seqline;
  }
}
close FULL;

#Parse OTU file and filter sequences
my %use_seq;
my $file_otu_filtered = $file_otu;
$file_otu_filtered =~ s/\.otu$/.filtered.otu/;
open (OLDOTU, $file_otu) or die;
open (NEWOTU, ">$file_otu_filtered") or die;
LINE:
while (my $line = <OLDOTU>) {
  chomp $line;
  if ($line =~ /^>/) {
    $line =~ s/^>//;
    my $curr_otu = "B_" . $line;
    my $print_otu = 0;
    my @curr_seq;
    while (my $otu_line = <OLDOTU>) {
      chomp $otu_line;
      if ($otu_line =~ /^$/) {
        (print NEWOTU ">$curr_otu\n" . (join "\n" => @curr_seq) . "\n\n") if $print_otu;
        next LINE;
      }
      if (exists $Seq{$otu_line}) {
        push @curr_seq => $otu_line;
        $use_seq{$otu_line} = 1;
        $print_otu = 1;
      } else {print "Removed sequence (not in alignment):\t$otu_line\n"}
    }
  }
}
close OLDOTU; close NEWOTU;

#Export filtered alignment
my $file_alignment_filtered = $file_alignment_full;
$file_alignment_filtered =~ s/\.fasta/.filtered.fasta/;
open(FILTERED, ">$file_alignment_filtered") or die;
foreach my $seq (nsort keys %use_seq) {
  print FILTERED ">$seq\n$Seq{$seq}\n"
}
