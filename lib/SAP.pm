#!/usr/local/bin/perl
##############################################################################
# PACKAGE NAME: SAP
# DESCRIPTION:  module for maniuplating protein sequence data
#               in the form of single amino acid polymorphism
#
# FUNCTIONS:
#   loadSAP - populate hash tables from .sap files
#       $locallele{$l}{$a}++;              A, 0101      -> A24
#       $locallelesap{$l}{$a}{$sap}++;     A, 0101, A24 -> 1
#       $locallelepos{$l}{$a}{$pos}=$sap;  A, 0101, 24  -> A24
#       $locsap{$l}{$sap}++;               A, A24       -> 1
#       $locpossap{$l}{$pos}{$sap}++;      A, 24, A24   -> 1
#
#      feat - make a feature vector for each allele 
#   seqDiff(loc, t1, t2) returns and array of pos:aa1:aa2 where alleles differ
#   humoralMisMatch(loc, d1, d2, r1, r2) returns hash of pos=>num_mismatches
#       based on truth table based on humoral (solid organ) rejection
#
#
# WRITTEN BY:   Martin Maiers
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ---------     ----------      -------------------------------------
# 1.000   2001-12-04     Maiers          Original code
# 1.001   2004-12-03     Maiers          Add properties
# 1.002   2010-12-03     Maiers          Add humoralMisMatch
#
##############################################################################
package SAP;
use strict;

#use vars qw (%locallele %locallelesap %locsap %locpossap %locallelepos %saphap);
use vars qw (%locallele %locallelesap %locsap %locpossap %locallelepos %struct  %Properties);
my @properties = ("charge", "size", "polarity");
my $init_struct = 0;
my $init_properties = 0;
my $hladb = "3.56.0";
my $exon3 = 0;
my $sapdir = "$ENV{'HOME'}/src/git/imgt2aa/data/";
my @loci = qw/HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1/;


sub set_exon3 {
  $exon3 = shift;
}
sub set_hladb {
  $hladb = shift;
}

##############################################################################
# function: feat
# parameters: loc, \@alleles, begin:end, begin:end, ...
# description: produce tables of features based on arguments
#   @alleles - restrict to positions polymorphic in these alleles
#              if empty then use all alleles
#   begin:end are positions, 
#   begin:end are positions, 
##############################################################################
sub feat {
  my ($loc, $ralleles, $rAF, @range) = @_;
  my $picminSAP = 0.0005;
  my $picminfeat = 0.0005;
  my $maxcombo = 1;

  my @ale = ();
  #
  # if a specific allele list is give, use that
  #
  if(@{$ralleles}) {
    @ale = @{$ralleles};
  } else {
    @ale = keys %{$locallelesap{$loc}};
  }

  
  my @p;  # polymorphic positions


  #
  # find the polymorphic positions to consider
  #
  foreach my $pos (keys %{$locpossap{$loc}}) {
    if (@range) {
      my $in =0;
      foreach my $r (@range) {
        my ($begin, $end) = split /:/, $r;
        $in++ if $pos >=$begin && $pos <= $end;
      }
      next unless $in;
    }
    push @p, $pos;
  }
  #
  # %saps contains all of the SAPs possible
  #
  my %saps;

  #
  # trim the list based on the PIC of the SAPs
  #
  if(@ale) {
    my @keep=();
    foreach my $pos (@p) {
      my %h;
      foreach my $allele (@{$ralleles}) {
        $h{$locallelepos{$loc}{$allele}{$pos}}++;
      }
      push @keep, $pos if PIC(\%h,$picminSAP) > $picminSAP;
    }
    @p = @keep;
  }
  
  
  #
  # %ps contains all positions and saps in those positions
  #
  my %ps;
  foreach my $pos (@p) {
    foreach my $allele (@ale) {
      my $sap = $locallelepos{$loc}{$allele}{$pos};
      $saps{$sap}++;
      $ps{$pos}{$sap}++;
    }
  }


  #
  # now that %saps is built, create a list of all possible
  # features
  #
  print scalar (@p), " positions: ", join (' ', sort {$a<=>$b} @p), "\n";
  print scalar (keys %saps), " saps: ", join (' ', keys %saps), "\n";
  my $nsaps = scalar keys %saps;
  
 
  #
  # now get all possible combinations of positions
  #
  my %combo;
  my @arr=selectN(\@p, $maxcombo);
  my %F;  # features	$F{"M45-S67"}{"0702"};
  foreach my $ch (@arr) {		# : delmited array of pos
    #
    # now get all possible combinations of saps at these  positions
    # ...and store the actuals
    #
    foreach my $allele (@ale) {
      my @forge=();			# its pronounced for-hay
      foreach my $pos (sort {$a<=>$b} split /:/, $ch) {		# a position
        my $sap = $locallelepos{$loc}{$allele}{$pos};
        push @forge, $sap if defined $sap;
      }
      $F{join('-', @forge)}{$allele}++;
    }
  }

  # %F{feature}{allele}
  #
  # prune the feature list based on polymorphism information content
  # take them in order of length; keep small features that hit the
  # same alleles as big ones
  #
  my %AAA;
  foreach my $f (sort {length($b)<=>length($a)} keys %F) {
    my %h;
    $h{1} = scalar (keys %{$F{$f}});	# the number of alleles with feature
    $h{0} = scalar (@ale) - $h{1};  	# the number of alleles without

    #
    # keep only those with unique allele sets
    #
    #my $al = join (':', keys %{$F{$f}});    
    #delete $F{$f}, next if defined $AAA{$al}; 
    #$AAA{$al}++;

    #
    # keep only those with sufficient infmation content
    #
    delete $F{$f}, next unless PIC(\%h,$picminfeat) > $picminfeat;
    #printf "%30s\t%s\n", $f, join (',', keys %{$F{$f}});
  }

 
  #
  # build the allele-feature output hash
  #
  my @features = (keys %F); 
  $$rAF{"ALL"} = join (':', @features);
  foreach my $allele (@ale) {
    my @af = ();
    foreach my $f (@features) {
      if (defined $F{$f}{$allele}) {
        push @af, "1";
      } else {  
        push @af, "0";
      }
    }
    $$rAF{$allele} = join (':', @af);
  }
  return;
}

##############################################################################
# function: chooseN
# parameters: $rarr ref to array of all positions
#             $rr   ref to array of combinations of positions : delim
#             $n    NNothing
##############################################################################
sub chooseN {
  my ($rarr, $rr, $nn) = @_;
  my $n = scalar (@{$rarr});	# stay away from the lion cage
  my $max = 2**$n;
  for(my $i=1; $i<$max; $i++) {
    my @in = ();
    my $b = numtobitstring($i);
    for(my $j=0; $j<$n; $j++) {
      push @in, $$rarr[$j] if (substr($b, $j, 1));
    }
    #push @{$rr}, join (':', @in) unless scalar(@in) >3;
    push @{$rr}, join (':', @in);
  }
}

sub selectN {
  my ($rarr, $n, $pos, $str) = @_;
  $pos = 0 unless defined $pos;
  
  #
  # recursion step
  #
  if($n<=0) {
    if(defined $str && length $str) {
      return "$str";
    } else {
      return;
    }
  }
  if($pos+1 > scalar(@{$rarr})) {
    if(defined $str && length $str) {
      return "$str";
    } else {
      return;
    }
  }


  my $newstr = "";
  my $item = $$rarr[$pos];
  if (!defined $str || !length $str) {
    $str = "";
    $newstr = $item;
  } else {
    $newstr="$str:$item";
  }
  
  # the  string
  # the recursive call if you add the current position
  # the recursive call if you don't add the current position
  my @recadd = selectN($rarr, $n-1, $pos+1, $newstr); 
  my @rec = selectN($rarr, $n, $pos+1, $str);
  if ($pos+1 > scalar(@{$rarr})) { 
    if (defined $str && length $str) {
      return ($str, @rec);
    } else {
      return (@rec);
    }
  } else {
    return (@rec, @recadd);
  }
}

##############################################################################
# function: chooseAll
##############################################################################
sub chooseAll {
  my ($rarr, $rr) = @_;
  my $n = scalar (@{$rarr});
  my $max = 2**$n;
  for(my $i=1; $i<$max; $i++) {
    my @in = ();
    my $b = numtobitstring($i);
    for(my $j=0; $j<$n; $j++) {
      push @in, $$rarr[$j] if (substr($b, $j, 1));
    }
    push @{$rr}, join (':', @in) unless scalar(@in) >3;
  }
}

sub numtobitstring {
  unpack("b*", pack("V", shift));  # pack in Vax order
}
##############################################################################
# function: PIC
##############################################################################
sub PIC {
  my $rh =shift;
  my $picmin =shift;
  $picmin = 0.01 unless defined $picmin; 
  my $sum=0;
  foreach (values %{$rh}) {
    $sum += $_;
  }
  my %P;
  foreach (keys %{$rh}) {
    $P{$_} = $$rh{$_}/$sum;
  }
  my $sum_pi_sq = 0;
  foreach my $i (keys %P) {
    $sum_pi_sq += $P{$i} * $P{$i};
  }
  my $sum_pi_pj_sq = 0;
  foreach my $i (keys %P) {
    foreach my $j (keys %P) {
      next if $i == $j;
      $sum_pi_pj_sq += $P{$i} * $P{$i} * $P{$j} * $P{$j};
    }
  }
  return (1 - $sum_pi_sq - $sum_pi_pj_sq);
}

##############################################################################
# function: loadSAP
##############################################################################
sub loadSAP {
  my $include_loc = shift;
  foreach my $loc (@loci) {
    my $file;
    if ($exon3) {
      $file = "$sapdir/$hladb.exon3/$loc.sap";
    } else {
      $file = "$sapdir/$hladb/$loc.sap";
    } 
    open FILE, $file or die "$0: $file $! \n";
    while(<FILE>) {
      chomp;
      my($l, $a, $s) = split /	/;
      my @saps = split / /, $s;
      foreach my $sap (@saps) {
        next unless $sap =~/(\D)(\d+)/;
        $sap = "$loc.$sap" if(defined $include_loc &&  $include_loc);
        my $pos = $2;
        $locallele{$l}{$a}++;
        $locallelesap{$l}{$a}{$sap}++;
        $locallelepos{$l}{$a}{$pos}=$sap;
        $locsap{$l}{$sap}++;
        $locpossap{$l}{$pos}{$sap}++;
#       foreach my $sap2 (@saps) {
#          $saphap{$sap}{$sap2}++;
#          $saphap{$sap2}{$sap}++;
#       }
      }
    }
  }
}
##############################################################################
# Function: isOnHap
# Description: return true if SAPS are on hap
##############################################################################
#sub isOnHap {
#  my ($sap1, $sap2)  = @_;
#  return 1 if defined $saphap{$sap1}{$sap2};
#  return 0;
#}

##############################################################################
# Function: hasSubFeat
# Description: return true if loc*type has subfeatures
##############################################################################
sub hasSubFeat {
  my ($allele, $subfeatures)  = @_;
  my ($loc, $typ) = split /\*/, $allele;
  my @sf = split /[\^:~]/, $subfeatures;
  my $ct=0;
  foreach my $subf (@sf) {
    $ct++ if $locallelesap{$loc}{$typ}{$subf};
  }  
  return 1 if $ct == scalar(@sf);
  return 0;
}

##############################################################################
# Function: hasSubFeatures
# Description: return true if type has subfeatures
##############################################################################
sub hasSubFeatures {
  my ($loc, $typ, $subfeatures)  = @_;
  my @sf = split /[\^:~]/, $subfeatures;
  my $ct=0;
  foreach my $subf (@sf) {
    $ct++ if $locallelesap{$loc}{$typ}{$subf};
  }  
  return 1 if $ct == scalar(@sf);
  return 0;
}

##############################################################################
# function: posStruct
##############################################################################
sub posStruct {
  my $pos = shift;
  loadStruct() unless $init_struct;
  return $struct{$pos} if defined $struct{$pos};
  warn "undefined position: $pos\n";
  return undef;
}

##############################################################################
# function: loadStruct
##############################################################################
sub loadStruct {
  $init_struct++;
  my $file = "$sapdir/class1_struct.txt";
  open FILE, $file or die "$0: $file $! \n";
  while(<FILE>) {
    chomp;
    my($category, $ranges) = split /:/;
    my $cat = substr($category, 0, 1);
    foreach (split /,/, $ranges) {
      my ($low, $high) = split /\-/;
      foreach ($low .. $high) {
        $struct{$_}=$cat;
      }
    }
  }
}

##############################################################################
# FUNCTION: seqDiff
# DESCRIPTION: output array of "44:R:K"
##############################################################################
sub seqDiff {
  my($loc, $typ1, $typ2) = @_;
  my @ret=();
  return undef unless defined $locallelepos{$loc}{$typ1};
  return undef unless defined $locallelepos{$loc}{$typ2};
  
  foreach my $pos (sort {$a<=>$b} keys %{$locallelepos{$loc}{$typ1}}) {
    my $sap1 = $locallelepos{$loc}{$typ1}{$pos};
    my $sap2 = $locallelepos{$loc}{$typ2}{$pos};
    
    next if $sap1 eq $sap2;
    #next if $sap1_s eq $sap2_s;
    #my ($loc1, $sap1) = split /\./, $sap1_s;
    #my ($loc2, $sap2) = split /\./, $sap2_s;
    my $aa1 = substr($sap1, 0, 1);
    my $aa2 = substr($sap2, 0, 1);

    push @ret, "$pos:$aa1:$aa2";
  }
  return @ret; 
}

##############################################################################
# FUNCTION: humoralMisMatch
# DESCRIPTION: loc, d1, d2, r1, r2 returns hash of pos=> num mismatches
##############################################################################
sub humoralMisMatch {
  my($loc, $d1, $d2, $r1, $r2, $rH) = @_;
  #print join ('|', "humoralMisMatch", $loc, $d1, $d2, $r1, $r2), "\n";

  #
  # homozygous 'hook'
  #
  my $d_hom = 0;
  if ($d1 eq $d2) {
    # donor homozygous; use truth table with max mismatch == 1
    $d_hom = 1;
  }

  foreach my $pos (sort {$a<=>$b} keys %{$locallelepos{$loc}{$d1}}) {
    my $d1aa = $locallelepos{$loc}{$d1}{$pos};
    my $d2aa = $locallelepos{$loc}{$d2}{$pos};
    my $r1aa = $locallelepos{$loc}{$r1}{$pos};
    my $r2aa = $locallelepos{$loc}{$r2}{$pos};

    my $m = getTruth($d1aa, $d2aa, $r1aa, $r2aa);
    $m=1 if ($m==2 && $d_hom ==1);
    $$rH{$pos} = $m;
  }
}


##############################################################################
# FUNCTION: getTruth
# DESCRIPTION: 2 donor and 2 recip amino acids; return 0, 1, 2 mismatches
##############################################################################
sub getTruth {
  my ($d1_aa, $d2_aa, $r1_aa, $r2_aa) = @_;
  my $m=0;
  $m++ if (defined $d1_aa && defined $r1_aa && $d1_aa ne $r1_aa) &&
    	  (defined $d1_aa && defined $r2_aa && $d1_aa ne $r2_aa);
  $m++ if (defined $d2_aa && defined $r1_aa && $d2_aa ne $r1_aa) &&
    	  (defined $d2_aa && defined $r2_aa && $d2_aa ne $r2_aa);
  #print join (' ', "GETTRUTH", $d1_aa, $d2_aa, $r1_aa, $r2_aa, $m), "\n";
  return $m;
}

##############################################################################
# FUNCTION: featDiffFirst
# DESCRIPTION: diff two alleles, output array of "B.S116, B.R82~B.G83""
#              which are the residues present in the first allele not the 2nd
##############################################################################
sub featDiffFirst {
  my($rpfeat,$loc, $typ1, $typ2) = @_;
  my ($l1, $t1) = split /\*/, $typ1;
  my ($l2, $t2) = split /\*/, $typ2;
  my @ret=seqDiff($loc, $t1, $t2);
  my @saps;
  foreach (@ret) {
    my ($pos, $res1, $res2) = split /:/;
    push @saps, "$loc.$res1"."$pos";
  }
  my %ret= ();
  foreach my $feat (@{$rpfeat}) {
    my @sf = split /[\^:~]/, $feat;
    foreach my $subfeature (@sf) {
      foreach my $sap (@saps) {
        $ret{$feat}++ if $subfeature eq $sap;
      }
    }
  }
  return sort keys %ret;
}

##############################################################################
# FUNCTION: featDiff
# DESCRIPTION: output array of "B.S116, B.R82~B.G83""
##############################################################################
sub featDiff {
  my($rpfeat,$loc, $typ1, $typ2) = @_;
  my ($l1, $t1) = split /\*/, $typ1;
  my ($l2, $t2) = split /\*/, $typ2;
  my @ret=seqDiff($loc, $t1, $t2);
  my @saps;
  foreach (@ret) {
    my ($pos, $res1, $res2) = split /:/;
    push @saps, "$loc.$res1"."$pos";
    push @saps, "$loc.$res2"."$pos";
  }
  my %ret= ();
  foreach my $feat (@{$rpfeat}) {
    my @sf = split /[\^:~]/, $feat;
    foreach my $subfeature (@sf) {
      foreach my $sap (@saps) {
        $ret{$feat}++ if $subfeature eq $sap;
      }
    }
  }
  return sort keys %ret;
}

##############################################################################
# function: loadProperties
##############################################################################
sub loadProperties {
  $init_properties++;
  my $namefile = "$sapdir/properties/name.txt";
  my %N;
  open FILE, $namefile or die "$0: $namefile $! \n";
  while(<FILE>) {
    chomp;
    my($name, $aa) = split /:/;
    $N{$name}= $aa;
  }

  foreach my $p (@properties) {
    my $file = "$sapdir/properties/$p.txt";
    open FILE, $file or die "$0: $file $! \n";
    while(<FILE>) {
      chomp;
      my($name, $value) = split /:/;
      die "unknown name: $name" unless defined $N{$name};
      $Properties{$p}{$N{$name}} = $value;
    }
    # normalize
    my $max_diff = 0;
    foreach my $aa1 (keys %{$Properties{$p}}) {
      foreach my $aa2 (keys %{$Properties{$p}}) {
        my $diff = abs ($Properties{$p}{$aa1} - $Properties{$p}{$aa2});
        $max_diff= $diff if $diff > $max_diff;
      }
    }
    foreach my $aa (keys %{$Properties{$p}}) {
      $Properties{$p}{$aa} = $max_diff ? $Properties{$p}{$aa}/$max_diff : 0;
    }
  }
  $init_properties = 1;
}
  

sub getProperty {
  loadProperties() unless $init_properties;
  my ($p, $aa) = @_;
  die "unknown property: $p" unless defined $Properties{$p};
  die "unknown AA: $aa" unless defined $Properties{$p}{$aa};
  return $Properties{$p}{$aa};
}

sub compareProperties {
  my ($p, $aa1, $aa2) = @_;
  my $p1 = getProperty($p, $aa1);
  my $p2 = getProperty($p, $aa2);
  return abs ($p1 - $p2);
}
1;
