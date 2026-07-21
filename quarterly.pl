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
use FindBin qw($Bin);
use Cwd qw(abs_path);

my $rootdir = abs_path($Bin);
my $datadir = "$rootdir/data";
my $hladb = "3.56.0";
my %options;
getopts("v:e", \%options);
$hladb = $options{v} if defined $options{v} && $options{v};

$hladb = "$hladb.exon3" if $options{e};
mkdir $datadir unless -e $datadir;
mkdir "$datadir/$hladb" unless -e "$datadir/$hladb";

chdir "$datadir/$hladb" or die "can't chdir to $datadir/$hladb: $!";
system("perl", "$rootdir/mksapTranslate") == 0 or die "mksapTranslate failed: $?";
system("perl", "$rootdir/mkkir") == 0 or die "mkkir failed: $?";

exit 0;
