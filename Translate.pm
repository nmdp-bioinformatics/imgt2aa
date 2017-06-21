#!/usr/local/bin/perl58
##############################################################################
# PACKAGE NAME:	Translate.pm
# DESCRIPTION:	AA translation from nucleotide sequence
# PARAMETERS:	none
# OUTPUT:	
# TABLES:       
#
# DATE WRITTEN: 2007-01-12
# WRITTEN BY:   Martin J. Maiers
#
# REVISION HISTORY: 
# REVISION DATE		REVISED BY	DESCRIPTION 
# ------- ----------	--------------	-------------------------------------
#
#       COPYRIGHT (C) 1994, 1995, 1997 NATIONAL MARROW DONOR PROGRAM.  
#                           ALL RIGHTS RESERVED        
##############################################################################
package Translate;
use strict;    # always
use warnings;  # or else


my $init=0;
my %tran;

sub initialize {
  $init=1;
  my $translate_dir = ".";
  my $prot="$translate_dir/protein.db";

  open(PROTEIN, $prot) or die "can't open $prot for reading";
  while(<PROTEIN>) {
    chomp;
    my ($aa,$aa_3,$aa_name,$na_list) = split /:/;
    my(@na) = split (/,/, $na_list);
    foreach my $n (@na) {
      $tran{$n}=$aa;
    }
  }
}

sub translate {
  my ($seq) = shift;
  my $ret = "";
  initialize() unless $init;
  my $length = length $seq;
  for(my $i=0; $i<$length; $i+=3) {
    my $codon = substr($seq,$i, 3);
    if (defined $tran{$codon}) {
      $ret.=$tran{$codon};
    } else {
      $ret.="-";
    }
  }
  return $ret;
} 
1;
