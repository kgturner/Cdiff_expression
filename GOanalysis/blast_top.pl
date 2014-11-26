################# tophit.pl ####################
# Created by Kay Hodgins
# 
#
# Assumes: blast file sorted with best hit at top of each list
# Returns: top hit for each contig
# Run: perl tophit.pl <blast file sorted> 
#	blast results are in tabluar format
#	
# 
#
########################################################
#!/usr/bin/perl
  use strict;
  use warnings;
#use lib '/Users/kayhodgins/lib';  #'/SciBorg/homes/hodgins/lib';#
#use Math::Random;



open(IN, "$ARGV[0]"); 
open(OUTFILE, ">unique.$ARGV[0]");


my @unique = ();
my %seen   = ();
my @split = ();

	foreach (<IN>){
		chomp $_;
		@split=split;
                next if $seen{$split[0]}++;
                print OUTFILE "$_\n";
	}

exit;




