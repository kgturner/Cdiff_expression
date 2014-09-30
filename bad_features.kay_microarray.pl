#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

my $xys = $ARGV[0]; #name of xys file
chomp $xys;
my $xmin = $ARGV[1]; #minimum x coordinate
chomp $xmin;
my $ymin = $ARGV[2]; #min y coord
chomp $ymin;
my $xmax = $ARGV[3]; #maximum x coordinate
chomp $xmax;
my $ymax = $ARGV[4]; #max y coord
chomp $ymax;

open(XYS, "<$xys") or die ("cannot open file: $xys\n"); #open xys file
open(OUTFILE, ">${ARGV[0]}_bfr"); #xys with bad features removed


#read in xys and replace values for positions in bad areas, and remove control probes
while (my $line = <XYS>) { 
    if ($line =~ /^#/) { #skip metadata line
	print OUTFILE $line;
	next;
    }
    my ($x, $y, $signal, $count, @arr) = split /\t/, $line; #split the line into an array at tabs
    if ($x eq 'X'){
	print OUTFILE $line;
	next;
    }
    my $positions = join(" ", $x,$y);
    if ($xmax >= $x && $x >= $xmin && $ymax >= $y && $y >= $ymin) { #do the replacement
	$signal = "NA";
	$count = "NA";
	print OUTFILE $x, "\t", $y, "\t", $signal, "\t", $count, "\n";
    }
    else { #this is an experimental probe in a 'safe' area of the array
	print OUTFILE $line;
    }
}
close XYS;
close OUTFILE;

exit;
