#!bin/perl

my $in=$ARGV[0];
my $out=$in . "_flip";
my %ortho=();

open IN, $in;
open OUT, ">", $out or die "Couldn't create out file";
while(<IN>){
		chomp;
		my @split=split;
		push @{$ortho{$split[0]}},$split[1];
}

foreach my $k (keys %ortho){
		print OUT "$k\t";
		my $len=scalar(@{$ortho{$k}});
		for my $g ( 0..($len-1)){
				if($g< $len-1){
						print OUT "$ortho{$k}[$g], ";
				}else{
						 print OUT "$ortho{$k}[$g]\n";
				}
		}
}

close IN;
close OUT;
