#!/usr/bin/perl
################################################################################
################################################################################
#Toolbox 16S: convenience script to translate MAPseq otutax files into OTU files.
##
#Usage:
#./translate.otutax_to_otu_files.pl <otu.tax>
#
#sebastian.schmidt@embl.de 2016-12-21
################################################################################
################################################################################

use strict;
use warnings;
use Sort::Naturally;

#Get input
my $file_otutax = shift @ARGV or die;

#Read otutax file and parse OTU mapping
my %OTU;
open(OTUTAX, $file_otutax) or die;
while (my $line = <OTUTAX>) {
  chomp $line;
  next if $line =~ /^#/;
  my ($seq, $otus_raw) = split /\t/ => $line;
  my ($domain, @otus) = split /;/ => $otus_raw;
  #Limit to bacteria (for now)
  next unless $domain eq "Bacteria";
  $OTU{F}{shift @otus}{$seq} = 1;
  $OTU{G}{shift @otus}{$seq} = 1;
  $OTU{S}{shift @otus}{$seq} = 1;
  $OTU{SS}{shift @otus}{$seq} = 1;
}
close OTUTAX;

#Export
foreach my $otu_lev ("F", "G", "S", "SS") {
  my $file_output = $file_otutax;
  $file_output =~ s/.otutax//;
  $file_output .= ".$otu_lev.otu";
  
  open (OTU, ">$file_output") or die;
  #Loop through OTUs and export
  foreach my $otu (nsort keys %{$OTU{$otu_lev}}) {
    print OTU ">$otu\n" . (join "\n" => (nsort keys %{$OTU{$otu_lev}{$otu}})) . "\n\n"
  }
  close OTU
}
