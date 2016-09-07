#!/usr/bin/perl
use strict;
use warnings;
use File::Glob qw(:globally :nocase);

my $filesDirectory = 'C:/Users/Kat/Documents/imputefeatures'; #Folder containing the xys files
my $knnFile = 'knn_replaced.txt'; #name of file with knn_replaced values
my $knnReplaceScript = 'knnreplace.kay_microarray.pl'; #name of script doing replacing
my $shortxysFile;
my @xysFileSplit;


#Run for each xys file in dorectory
while (my $xysFile = <$filesDirectory/*.xys>)  
{
	@xysFileSplit = split(/\//, $xysFile);
	$shortxysFile = $xysFileSplit[-1];
	#print "Command will be:\nperl $knnReplaceScript $shortxysFile $knnFile\n";
	#die;
	`perl $knnReplaceScript $shortxysFile $knnFile`;
}

