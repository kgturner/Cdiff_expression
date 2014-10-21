#!/usr/bin/perl
use strict;
use warnings;

open INFILE, "<", $ARGV[0] or die $!;

while(my $line = <INFILE>)
{
	chomp($line);
	`perl bad_features.kay_microarray.pl $line`;
}
