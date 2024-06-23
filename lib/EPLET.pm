#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:	EPLET.pm
# DESCRIPTION:	load eplet registry
#
# DATE WRITTEN: 2024-06-18
# WRITTEN BY:   Martin Maiers
#
############################################################################
package EPLET;
use strict;    # always
use warnings;  # or else
use Config::JSON;
use lib "$ENV{'HOME'}/src/git/imgt2aa/lib";
use SAP;


my $init = 0;
my %E;


sub initialize  {
  my $config_file  = shift;
  my $config = Config::JSON->new($config_file);

  SAP::loadSAP();
  my $epregistry= $config->get("epregistry");
  my @epsets = qw/ABC DP DQ DRB/;

  foreach my $epset (@epsets) {
    my $file = "$epregistry/HLA-$epset.csv";
    open FILE, $file or die "$!: $file";
    while (<FILE>) {
      #ID,Name,Description,Exposition ,3D View,Confirmation,Evidence ,Str. Epitope,Frequency,Luminex? Alleles,All Alleles
      chomp;
      my ($id, $name, $description, $exposition, $view, $confirmation, $dvidence, $str, $freq, $lum_alleles, $all_alleles) =split /,/;
      next if $description =~/Description/;
      $E{$epset}{$id}{name} = $name;
      $E{$epset}{$id}{description} = $description;
      $E{$epset}{$id}{subfeat} = join (':', getsubfeat($description));
      $E{$epset}{$id}{exposition} = $exposition;
      $E{$epset}{$id}{all_alleles} = $all_alleles;
    }
  }
  $init++;
}

sub getsubfeat {
  my $s = shift;
  my @ret = ();
  while ($s =~/(\d+)(\D+)/g) {
    my $startpos = $1; 
    my $aas = $2; 
    # remove + and - because they have no meaning
    $aas=~s/[\+\-]//g;
    
    # offset is for multi-position eplets (e.g. 65QIT)
    my $offset = 0;
    while ($aas =~/(\D)/g) {
      my $aa = $1; 
      next unless $aa=~/[A-Z]/;
      my $pos = $startpos + $offset;
      #print "pos: $pos  AA: $aa\n";
      push @ret, $aa.$pos;
      $offset++;
    }   
  }
  return @ret;
}
sub get {
  die "initialize first" unless $init;

  # put allele into WHO format
  my $allele = get_who(shift);
  my $loc = (split /\*/, $allele)[0];
  my $epset = get_epset($loc);
  if (!defined $epset) {
    die "no eplets for loc $loc";
  }


  my @ret=();
  foreach my $id (keys %{$E{$epset}}) {
    if (SAP::hasSubFeat($allele, $E{$epset}{$id}{"subfeat"})) {
      push @ret, $E{$epset}{$id}{"name"};
    }
  }
  return sort epsort @ret;
}


sub get_epset {
  my $loc = shift;
  return "ABC" if $loc =~/^[ABC]/;
  return "ABC" if $loc =~/^HLA-[ABC]/;
  return "DQ" if $loc =~/^DQ/;
  return "DQ" if $loc =~/^HLA-DQ/;
  return "DP" if $loc =~/^DP/;
  return "DP" if $loc =~/^HLA-DP/;
  return "DRB" if $loc =~/^DRB/;
  return "DRB" if $loc =~/^HLA-DRB/;
  return undef;
}

sub mismatches {
  my($rd, $rr) = @_;
  my %disjoint;
  my %D;
  my %R;
  foreach my $eplet (@{$rr}) {
    $R{$eplet}++;
  }
  foreach my $eplet (@{$rd}) {
    $D{$eplet}++;
  }

  # donor (HvG) only
  foreach my $eplet (keys %D) {
    $disjoint{$eplet}++ unless defined $R{$eplet};
  }
  #foreach my $eplet (keys %R) {
  #  $disjoint{$eplet}++ unless defined $D{$eplet};
  #}
  return (sort epsort keys %disjoint);
}

sub epsort {
  $a=~/(\d+)(\D+)/;
  my $da = $1;
  my $Da = $2;
  $b=~/(\d+)(\D+)/;
  my $db = $1;
  my $Db = $2;
  if ($da < $db) {
    return -1;
  } elsif($da == $db) {
    return $Da cmp $Db;
  } else {
    return 1;
  }
}
sub get_who {
  my $s = shift;
  $s = "HLA-".$s unless $s=~/^HLA-/;
  return $s;
}
1;

