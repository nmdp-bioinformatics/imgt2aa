#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	quarterly.pl
# DESCRIPTION:	run the quarterly files
# PARAMETERS:	none
# OUTPUT:	
# TABLES:       
#
# DATE WRITTEN: 2020-08-31
# WRITTEN BY:   Martin Maiers
#
# REVISION HISTORY: 
# REVISION DATE		REVISED BY	DESCRIPTION 
# ------- ----------	--------------	-------------------------------------
#
#       COPYRIGHT (C) 2020 NATIONAL MARROW DONOR PROGRAM.  
#               ALL RIGHTS RESERVED        
##############################################################################
use strict;    # always
use warnings;  # or else
use Getopt::Std;

my $datadir = "data";
my $hladb = "3.56.0";
my %options;
getopts("v:e", \%options);
$hladb = $options{v} if defined $options{v} && $options{v};

$hladb = "$hladb.exon3" if $options{e};
mkdir $datadir unless -e $datadir;
mkdir "$datadir/$hladb" unless -e "$datadir/$hladb";

`cd "$datadir/$hladb"; perl ../../mksapTranslate ; perl ../../mkkir`;

exit 0;

