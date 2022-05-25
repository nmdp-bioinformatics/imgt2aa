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
use XML::Simple;
use lib ".";
use Translate;  # every Bioinformatician has written one

my $verbose =1;

# 
# download latest hla.xml 
# $ curl ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip >hla.xml.zip
# $ unzip hla.xml.zip
#
my $file = `cat hla.xml`;

print STDERR "parsing..." if $verbose;
my $alleles = XMLin($file);
print STDERR "done\n" if $verbose;

my %GF;  # Gene Feature

print STDERR "number of alleles: ", scalar keys %{$alleles->{"allele"}}, "\n";
foreach my $allele (keys %{$alleles->{"allele"}}) {
  # $allele='HLA-A*01:01:08'
  my $a = $alleles->{allele}->{$allele};
  my $s = $a->{sequence};
  my $nucsequence= $s->{nucsequence};
  foreach my $feature (keys %{$s->{feature}}) {
    next unless $feature =~/(Exon|Intron|UTR)/;
    # gDNA
    my $start = $s->{feature}->{$feature}->{SequenceCoordinates}->{start};
    my $end = $s->{feature}->{$feature}->{SequenceCoordinates}->{end};
    # cDNA
    my $rf = $s->{feature}->{$feature}->{cDNACoordinates}->{readingframe};
    my $cdna_start = $s->{feature}->{$feature}->{cDNACoordinates}->{start};
    my $cdna_end = $s->{feature}->{$feature}->{cDNACoordinates}->{end};
    # sequence of this feature
    my $nseq = substr($nucsequence, $start-1, $end-$start+1);
    $GF{$allele}{$feature}{nuc} = $nseq;
    $GF{$allele}{$feature}{rf}  = $rf;
    $GF{$allele}{$feature}{cdna_start}  = $cdna_start;
  }
}

#
# cDNA position corresponding to codon 1 in the mature protein
# 
my %PM;
$PM{"HLA-DPB1"} = 84;
$PM{"HLA-DQB1"} = 84;
$PM{"HLA-DRB1"} = 84;

# output hash
my %O; 

foreach my $allele (sort keys %GF) {
  my $loc = (split /\*/, $allele)[0];
  next unless defined $PM{$loc};

  foreach my $feature (sort keys %{$GF{$allele}}) {
    if ($loc=~/HLA-D/) {
      next unless $feature=~/Exon 2/;
    } else {
      next unless $feature=~/Exon [23]/;
    }
    my $rf = $GF{$allele}{$feature}{rf};
    my $nuc = $GF{$allele}{$feature}{nuc};
    my $aaseq = TranslateRF($nuc, $rf);
    #
    # figure out how much to left pad with asterisks based on cdna_start
    #
    my $cdna_start = $GF{$allele}{$feature}{cdna_start};
    my $pos1_mature= $PM{$loc};

    # this line is tricky
    # from cDNA start, substract the cDNA position of the codon 1 
    # then add rf this is the postion of the first codon, in frame
    # subtract 1 to get the offset (number of unspecified AA)
    my $offset = int ((($cdna_start+$rf)-$pos1_mature)/3)-1;
    print STDERR "$loc offset $offset\n";

    # fill in with this many asterisks
    my $prefix = "*" x $offset;
    # and prefix
    $aaseq = $prefix . $aaseq;

    # get the short name: DPB1*02:01 (from DPB1*02:01:01:01)
    my $allele2 = getallele2($allele);

    # store the result
    $O{$loc}{$allele2} = $aaseq;
  }
}

#
# print out
#
foreach my $loc (sort keys %O) {
  open LOCDB, ">$loc.db" or die "$!: $loc.db";
  foreach my $allele (sort keys %{$O{$loc}}) {
    my ($loc, $a) = split /\*/, $allele;
    my $shortloc = (split /\-/, $loc)[1];
    print LOCDB join ('	', $shortloc, $a, $O{$loc}{$allele}), "\n";
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
