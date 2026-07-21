#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	imgt2aa.pl
# DESCRIPTION:	parse the IMGT HLA.xml and generate a table of loc, allele, AA on stdout
#
# DATE WRITTEN: 2017-06-20
# WRITTEN BY:   Martin Maiers
##############################################################################
use strict;    # always
use warnings;  # or else
use XML::Twig;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use lib $Bin;
use Translate;  # every Bioinformatician has written one
use Getopt::Std;

my $rootdir = abs_path($Bin);
my $verbose =1;
my $hladb="3.65.0";
my %options;
getopts("ev:", \%options);
$hladb = $options{v} if defined $options{v} && $options{v};
my $datadir="$rootdir/data/$hladb";
`mkdir -p $datadir`;

my $exon3 = $options{e};
# 
# download latest hla.xml 
# $ curl ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip >hla.xml.zip
# $ unzip hla.xml.zip
#
my %GF;  # Gene Feature
my $allele_count = 0;
my $progress_interval = 1000;

print STDERR "parsing..." if $verbose;
my $twig = XML::Twig->new(
  twig_handlers => {
    allele => sub {
      my ($twig, $allele_elt) = @_;
      my $allele = $allele_elt->att('name');
      my $sequence_elt = $allele_elt->first_child('sequence');
      if (!$allele || !$sequence_elt) {
        $allele_elt->purge;
        return;
      }
      my $nucsequence = $sequence_elt->first_child_text('nucsequence') // "";
      if (!$nucsequence) {
        $allele_elt->purge;
        return;
      }

      $allele_count++;
      if ($verbose && $allele_count % $progress_interval == 0) {
        print STDERR "\rparsed $allele_count alleles...";
      }
      foreach my $feature_elt ($sequence_elt->children('feature')) {
        my $feature = $feature_elt->att('name') // "";
        next unless $feature =~ /(Exon|Intron|UTR)/;

        my $seq_coords = $feature_elt->first_child('SequenceCoordinates');
        next unless $seq_coords;
        my $start = $seq_coords->att('start');
        my $end = $seq_coords->att('end');
        next unless defined $start && defined $end;

        my $cdna_coords = $feature_elt->first_child('cDNACoordinates');
        my $rf = defined $cdna_coords ? $cdna_coords->att('readingframe') : undef;
        my $cdna_start = defined $cdna_coords ? $cdna_coords->att('start') : undef;

        my $nseq = substr($nucsequence, $start - 1, $end - $start + 1);
        $GF{$allele}{$feature}{nuc} = $nseq;
        $GF{$allele}{$feature}{rf} = $rf if defined $rf;
        $GF{$allele}{$feature}{cdna_start} = $cdna_start if defined $cdna_start;
      }

      $allele_elt->purge;
    },
  },
);
$twig->parsefile("$rootdir/hla.xml");
print STDERR "\rparsed $allele_count alleles...done\n" if $verbose;
print STDERR "number of alleles: $allele_count\n";

#
# cDNA position corresponding to codon 1 in the mature protein
# 
my %PM;
$PM{"HLA-DPB1"} = 84;
$PM{"HLA-DQB1"} = 93;
$PM{"HLA-DRB1"} = 84;
$PM{"HLA-A"} = 69;
$PM{"HLA-B"} = 69;
$PM{"HLA-C"} = 69;

# output hash
my %O; 

foreach my $allele (sort keys %GF) {
  my $loc = (split /\*/, $allele)[0];
  next unless defined $PM{$loc};

  my $nuc = "";
  foreach my $feature (sort keys %{$GF{$allele}}) {
    if ($loc=~/HLA-D/) {
      if ($exon3) {
        next unless $feature=~/Exon [1234]/;
      } else {
        next unless $feature=~/Exon [2]/;
      }
    } else {
      if ($exon3) {
        next unless $feature=~/Exon [234]/;
      } else {
        next unless $feature=~/Exon [23]/;
      }
    }
    if ($feature eq "Exon 1" && $exon3) {
      if ($loc=~/HLA-D[RQ]B[1345]/) {  
        $nuc .= substr($GF{$allele}{$feature}{nuc}, -13);
      }
    } else {
      $nuc .= $GF{$allele}{$feature}{nuc};
    }
  }

  my $rf = $GF{$allele}{"Exon 2"}{rf};

  # readingframe = "3" means this G|AG

  my $aaseq;
  if ($loc=~/HLA-[ABC]/) {
    my $exon1 = defined $GF{$allele}{"Exon 1"}{nuc} ?  $GF{$allele}{"Exon 1"}{nuc} : "*";
    my $exon4 = defined $GF{$allele}{"Exon 4"}{nuc} ?  $GF{$allele}{"Exon 4"}{nuc} : "**";
    $nuc= substr($exon1, -1, 1) . $nuc . substr($exon4, 0, 2);
    $aaseq = TranslateRF($nuc, 1);
  } else {
    # class II
    if ($exon3) {
      if (defined $GF{$allele}{"Exon 1"}{nuc}) {
        $aaseq = TranslateRF($nuc, 1);
      } else {
        $aaseq = TranslateRF($nuc, $rf);
      }
    } else {
      $aaseq = TranslateRF($nuc, $rf);
    }
  }


  if (!$exon3 || ! defined $GF{$allele}{"Exon 1"}{nuc}) { 
    #
    # figure out how much to left pad with asterisks based on cdna_start
    #
    my $cdna_start = $GF{$allele}{"Exon 2"}{cdna_start};
    my $pos1_mature= $PM{$loc};

    # this line is tricky
    # from cDNA start, substract the cDNA position of the codon 1 
    # then add rf this is the postion of the first codon, in frame
    # subtract 1 to get the offset (number of unspecified AA)
    my $offset = int ((($cdna_start+$rf)-$pos1_mature)/3)-1;

    if ($loc=~/HLA-[ABC]/) {
      $offset = int ((($cdna_start)-$pos1_mature)/3)-1;
    }
    if ($offset < 0) {
      print STDERR "$allele has negative offset $offset; clamping to 0\n" if $verbose > 1;
      $offset = 0;
    }
    print STDERR "$loc offset $offset\n" if $verbose > 1;
  
    # fill in with this many asterisks
    my $prefix = "*" x $offset;
  
    $aaseq = $prefix . $aaseq;
  } 

  # get the short name: DPB1*02:01 (from DPB1*02:01:01:01)
  my $allele2 = getallele2($allele);
  # store the result
  $O{$loc}{$allele2} = $aaseq;
}

#
# print out
#
foreach my $loc (sort keys %O) {
  if ($exon3) {
    open LOCDB, ">$datadir.exon3/$loc.db" or die "$!: $loc.db";
  } else {
    open LOCDB, ">$datadir/$loc.db" or die "$!: $loc.db";
  }
  foreach my $allele (sort keys %{$O{$loc}}) {
    my ($loc, $a) = split /\*/, $allele;
    my $shortloc = (split /\-/, $loc)[1];
    print LOCDB join ('	', $loc, $a, $O{$loc}{$allele}), "\n";
  }
}

exit 0;


# 
# get the 2-field allele name
#

sub getallele2 {
  my $la = shift;
  my ($loc, $allele) = split /\*/, $la;
  my (@af) = split /\:/, $allele;
  my $allele2 = join (':', $af[0], $af[1]);
  return join ('*', $loc, $allele2);
}

# 
# translate based on reading frame offset
#
sub TranslateRF {
  my ($nuc, $rf) = @_;
  if (!defined $rf || !$rf || $rf == 1) {
    return Translate::translate($nuc);
  } elsif ($rf ==2) {
    return Translate::translate(substr($nuc,1));
  } elsif ($rf ==3) {
    return defined substr($nuc,2) ? Translate::translate(substr($nuc,2)) : "";
  }
}
__END__
