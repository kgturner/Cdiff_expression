#!/usr/bin/perl
use strict;
use warnings;

#this script takes in an xys file for a nimblegen microarray, and a matrix containing imputed values for the same microarray in one of the columns
#the output is a new xys file with imputed values for NA spots that satisfied the criteria for the imputation method used

my $xys = $ARGV[0]; #name of 'corrected' xys file
chomp $xys;
my $knn = $ARGV[1]; #name of file with knn_replaced values
chomp $knn;

open(XYS, "$xys") or die ("cannot open $xys\n"); #xys file
open(KNN, "$knn") or die ("cannot open $knn\n"); #knn_replaced file
open(OUTFILE, ">${xys}_knn"); #xys with bad features removed

my $a=0;
my $j;
my %knndata;
while (my $line = <KNN>) {
    my @arr = split /\t/, $line;
    my $n = 0;
    my $length = @arr;
    if ($a == 0) {#figure out which element in the header of the KNN_REPLACE file matches the XYS file we are replacing the values of
	foreach my $b (@arr) {
	    chomp $b;
	    if (index($xys, $b)!= -1) { #if the name of the xys file and the column header match, that's the column number we want
		$j = $n;
		print STDERR "$b is column $j in the file\n";
	    }
	    $n++; #count of column headers checked so far
	}
	$a++; #count out of this loop. don't need to count all the lines though, just the header.
    }
    $knndata{$arr[0]}=$arr[$j]; #hash KNN_REPLACE FILE USING PROBE IDs
    #print STDERR "$arr[0] is the probe and $arr[$j] is the value\n";
}
close KNN;

#READ IN XYS FILE. COUNT LINE NUMBERS. REPLACE EVERY VALUE WITH AN N-2 VALUE FROM KNN_REPLACE FILE, print out lines after replacing values. done.
$a=0; #need to count lines again
while (my $line = <XYS>) {
    if ($a < 2) { #print out first two lines (header lines) without changing them
	print OUTFILE $line;
    }
    else {
	my ($x, $y, $val, $thing) = split /\t/, $line;
	my $c= $a-1;
	chomp $val;
	if ($val eq "NA") {
	    #print "$val\n";
	    if (exists $knndata{$c}){
		print STDERR "replacing $val with $knndata{$c}\n";
		$val=$knndata{$c};
		$thing=1;
	    }
	}
	print OUTFILE $x, "\t", $y, "\t", $val, "\t", $thing;
    }
    $a++;
}
close XYS;
close OUTFILE;
