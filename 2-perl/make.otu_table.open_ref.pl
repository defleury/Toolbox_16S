#!/usr/bin/perl
################################################################################
################################################################################
#Toolbox 16S: convenience script to process MAPseq/hpc-clust output and generate
#an open-reference OTU table.
#
#This script generates an open-reference OTU table and corresponding OTU
#taxonomic data from MAPseq output and an hpc-clust OTU file.
#
#Usage:
#./make.otu_table.open_ref.pl
#  <sample.mapping>              #input: sequence-to-sample mapping file
#  <MAPseq.ncbitax>              #input: MAPseq taxonomy mapping file
#  <MAPseq.otutax>               #input: MAPseq OTU mapping file
#  <denovo.OTUs>                 #input: de novo OTU clustering file
#  <level>                       #parameter: MAPseq OTU mapping level (F, G, R, S, SS)
#  <tax.consensus>               #parameter: threshold for consensus taxonomy
#  <min.otu_size>                #parameter: minimum OTU size (globally)
#  <OTU.table>                   #output: combined OTU table
#  <OTU.data>                    #output: combined OTU consensus taxonomies
#
#sebastian.schmidt@embl.de 2016-12-15
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
my $file_MAPseq_mapping_otutax = shift @ARGV or die;
my $file_OTU_denovo = shift @ARGV or die;
my $level = shift @ARGV or die;
my $thresh_tax_consensus = shift @ARGV or die;
my $thresh_min_otu_size = shift @ARGV or die;
my $file_OTU_table = shift @ARGV or die;
my $file_OTU_data = shift @ARGV or die;
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
my @tax_levels = ("phylum", "class", "order", "family", "genus", "species");

#Preallocate data structures
my %Sequence_per_Sample;
my %Sample_per_Sequence;
my %seq_tax;
my %OTU_Mapping;
my %OTU_data;

#Get parameters
my $lev_idx;
if ($level eq "F") {$lev_idx = 0}
elsif ($level eq "G") {$lev_idx = 1}
elsif ($level eq "S") {$lev_idx = 2}
elsif ($level eq "SS") {$lev_idx = 3}
else {die "ERROR: Unknown OTU clustering level: $level\n"}
##############################

##############################
#Read sequence to sample mapping
print "Reading sequence-to-sample mapping...";
open(MAP, $file_sample_mapping) or die;
while (<MAP>) {chomp; my ($seq, $smpl) = split /\t/; $Sequence_per_Sample{$seq} = $smpl; $Sample_per_Sequence{$smpl}{$seq} = 1}
close MAP;
my @ordered_samples = nsort keys %Sample_per_Sequence;
print "done.\n";
##############################

##############################
#Read sequence taxonomy
print "Reading sequence taxonomy...";
open (SEQTAX, $file_MAPseq_mapping_ncbitax) or die;
while (<SEQTAX>) {
   chomp; my ($acc, @tax, $consensus_tax) = split /\t/;
   foreach my $lev (@tax_levels) {
      my $tax = shift @tax;
      last if $tax ~~ ["", " ", "NA"];
      #Fix taxonomic artifact "PHY_Coriobacteriia"
      $tax =~ s/PHY_Coriobacteriia/Actinobacteria/;
      $seq_tax{$acc}{$lev} = $tax
   }
}
close SEQTAX;
print "done.\n";
##############################

##############################
#Read current reference OTU mapping
#=> filter by taxonomy (remove archaea, eukarya and unclassified at domain level)
##############################
print "Reading reference OTU mapping...";
open(OTUMAP, "$file_MAPseq_mapping_otutax") or die;
LINE:
while (<OTUMAP>) {
   chomp; my ($acc, @mapping) = split /\t/;
   my $tmp_otu = $mapping[$lev_idx];
   
   #Generate proper OTU name, by pre-fixing with domain and level
   my $curr_otu = "B_$level" . "_$tmp_otu";
   
   #Skip if current sequence did not map at current level
   next LINE if $curr_otu ~~ ["", " "];
   
   #Skip if current sequence was not mapped to a "known" OTU
   next LINE if $curr_otu eq "NA";
   
   #Get sample for current sequence
   my $sample = $Sequence_per_Sample{$acc};
   
   #Store current OTU
   $OTU_Mapping{$curr_otu}{Samples}{$sample}++;
   $OTU_data{$curr_otu}{Size}++;
   #Collect current taxonomy
   LEVEL:
   foreach my $tax_lev (@tax_levels) {
      last LEVEL unless exists $seq_tax{$acc}{$tax_lev};
      $OTU_Mapping{$curr_otu}{Collect_Taxonomy}{$tax_lev}{$seq_tax{$acc}{$tax_lev}}++;
   }
}
close OTUMAP;
print "done.\n";
##############################
#Read and process current de novo OTUs
##############################
print "Reading de novo OTU mapping...";
open (DNMAP, "$file_OTU_denovo") or die;
LINE:
while (my $otu_name = <DNMAP>) {
   if ($otu_name =~ /^>/) {
      chomp $otu_name; $otu_name =~ s/^>//;
      
      #Preallocate size counter
      my $size = 0;
      
      #Iterate through sequences assigned to current OTU
      SEQLINE:
      while (my $acc = <DNMAP>) {
         chomp $acc;
         
         #Process and store current OTU data when empty line is encountered
         if ($acc =~ /^$/) {
            next LINE unless $size > 0;
            $OTU_data{$otu_name}{Size} = $size;
            next LINE;
         }
         
         #Get sample for current sequence
         my $sample = $Sequence_per_Sample{$acc};
         #Store sample mapping
         $OTU_Mapping{$otu_name}{Samples}{$sample}++;
         $size++;
         
         #Collect current taxonomy
         LEVEL:
         foreach my $tax_lev (@tax_levels) {
            last LEVEL unless exists $seq_tax{$acc}{$tax_lev};
            $OTU_Mapping{$otu_name}{Collect_Taxonomy}{$tax_lev}{$seq_tax{$acc}{$tax_lev}}++;
         }
      }
   }
}
close DNMAP;
print "done.\n";
##############################

##############################
#Loop through OTUs to establish consensus taxonomy
print "Determining OTU consensus taxonomy...";
OTU:
foreach my $otu (keys %OTU_Mapping) {
   #Get current taxonomy collector hash and OTU size
   my %collect_tax = %{$OTU_Mapping{$otu}{Collect_Taxonomy}};
   my $size = $OTU_data{$otu}{Size};
   
   next OTU unless $size > 0;
   
   my $tax_string = "";
   
   #Get consensus taxonomy
   LEVEL:
   foreach my $tax_lev (@tax_levels) {
      if (exists $collect_tax{$tax_lev}) {
         my $largest = largest_value %{$collect_tax{$tax_lev}};
         my $dominance = $collect_tax{$tax_lev}{$largest} / $size;
         
         if ($dominance > $thresh_tax_consensus) {
            $OTU_data{$otu}{Consensus_Taxonomy}{$tax_lev} = $largest;
            $tax_string .= ";$largest";
         }
         else {last LEVEL}
      }
   }
   
   $tax_string =~ s/^;//;
   
   $OTU_data{$otu}{Consensus_Taxonomy} = $tax_string;
}
print "done.\n";
##############################

##############################
#Export
#=> OTU table
#=> OTU data
##############################
open (TABLE, ">$file_OTU_table") or die;
open (DATA, ">$file_OTU_data") or die;
#Preallocate OTU table (header)
print TABLE "\t" . join ("\t", @ordered_samples) . "\n";
#Print OTU data header
print DATA "OTU_name\tsize\tconsensus_taxonomy\n";
#Loop through OTUs and export data
print "Writing OTU table, OTU data and OTU file...";
foreach my $otu (nsort keys %OTU_Mapping) {
   #Skip current OTU if too few global reads attracted
   next unless $OTU_data{$otu}{Size} >= $thresh_min_otu_size;
   
   #Export to OTU table
   print TABLE "$otu";
   foreach my $smpl (@ordered_samples) {
      my $count = ($OTU_Mapping{$otu}{Samples}{$smpl} || 0);
      print TABLE "\t$count"
   }
   print TABLE "\n";
   
   #Export OTU data (consensus taxonomy)
   print DATA "$otu\t$OTU_data{$otu}{Size}\t$OTU_data{$otu}{Consensus_Taxonomy}\n";
}
close TABLE; close DATA;
print "done.\n";
##############################

