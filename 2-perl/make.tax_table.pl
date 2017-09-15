#!/usr/bin/perl
################################################################################
################################################################################
#Toolbox 16S: convenience script to generate taxonomy tables from MAPseq output
#
#This script parses a MAPseq taxonomy output file and generates taxon-to-sample
#mapping tables per taxonomic level
#
#Usage: make.tax_table.pl
#  <sample.mapping>              #input: sequence-to-sample mapping file
#  <MAPseq.ncbitax>              #input: MAPseq taxonomy mapping file
#  <thresh.tax_size>             #input: minimum number of seqs per taxon
#  <out.dir>                     #output: target directory to store tax tables
#
#sebastian.schmidt@embl.de 2017-09-15
################################################################################
################################################################################

use strict;
use warnings;
use Sort::Naturally;
$|++;

##############################
#Get input
my $file_sample_mapping = shift @ARGV or die;
my $file_MAPseq_mapping_ncbitax = shift @ARGV or die;
my $thresh_min_tax_size = shift @ARGV or die;
my $folder_output = shift @ARGV or die;
##############################

##############################
#Get a fast subroutine for key lookup by max value
sub largest_value (\%) {
      my $hash   = shift;
      my ($key, @keys) = keys   %$hash;
      my ($big, @vals) = values %$hash;

      for (0 .. $#keys) {
            if ($vals[$_] > $big) {
                  $big = $vals[$_];
                  $key = $keys[$_];
            }
      }
      $key
}
##############################

##############################
#Preallocate taxonomic levels
my @tax_levels = ("kingdom", "phylum", "class", "order", "family", "genus", "species");

#Preallocate data structures
my %Sequence_per_Sample;
my %Sample_Tax;
my %Taxa;
my %Samples;
##############################

##############################
#Read sequence to sample mapping
print "Reading sequence-to-sample mapping...";
open(MAP, $file_sample_mapping) or die;
while (<MAP>) {chomp; my ($seq, $smpl) = split /\t/; $Sequence_per_Sample{$seq} = $smpl; $Samples{$smpl}=1}
close MAP;
my @ordered_samples = nsort keys %Samples;
print "done.\n";
##############################

##############################
#Read sequence taxonomy
print "Reading sequence taxonomy...";
open (SEQTAX, $file_MAPseq_mapping_ncbitax) or die;
while (<SEQTAX>) {
   chomp; my ($acc, @tax, $consensus_tax) = split /\t/;
   #Get current sample
   my $smpl = $Sequence_per_Sample{$acc};
   foreach my $lev (@tax_levels) {
      my $tax = shift @tax;
      #Fix taxonomic artifact "PHY_Coriobacteriia"
      #=> these are actually Actinobacteria...
      $tax =~ s/PHY_Coriobacteriia/Actinobacteria/;
      #Relabel "NA" as "Unclassified"
      if ($tax eq 'NA') {$Sample_Tax{$lev}{$smpl}{Unclassified}++}
      else {$Sample_Tax{$lev}{$smpl}{$tax}++; $Taxa{$lev}{$tax}++}
   }
}
close SEQTAX;
print "done.\n";
##############################

##############################
#Export
#=> tax tables per taxonomic level
##############################
foreach my $lev (@tax_levels) {
   #Skip trivial level "kingdom"
   next if $lev eq "kingdom";
   #Set current tax table file name
   my $file_table = "$folder_output/tax_table.$lev.tsv";
   #Get taxa for current level
   my @ordered_taxa = ("Unclassified", nsort keys %{$Taxa{$lev}});
   $Taxa{$lev}{Unclassified} = $thresh_min_tax_size + 1;
   #Open output file
   open (TABLE, ">$file_table") or die;
   #Preallocate tax table (header)
   print TABLE "\t" . join ("\t", @ordered_samples) . "\n";
   #Loop through taxa and export data
   print "Writing $lev tax table...";
      #Skip current OTU if too few global reads attracted
   foreach my $tax (@ordered_taxa) {
      next unless $Taxa{$lev}{$tax} >= $thresh_min_tax_size;
      
      #Export to taxa table
      print TABLE "$tax";
      foreach my $smpl (@ordered_samples) {
         my $count = ($Sample_Tax{$lev}{$smpl}{$tax} || 0);
         print TABLE "\t$count"
      }
      print TABLE "\n";
   }
   close TABLE;
   print "done.\n";
}
##############################
