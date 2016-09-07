#!bin/perl

my $in=$ARGV[0];

my %ortho=();

open IN, $in;
while(<IN>){
		chomp;
		my @split=split;
		push @{$ortho{$split[0]}},$split[1];
}

foreach my $k (keys %ortho){
		print "$k\t";
		my $len=scalar(@{$ortho{$k}});
		for my $g ( 0..($len-1)){
				if($g< $len-1){
						print "$ortho{$k}[$g], ";
				}else{
						 print "$ortho{$k}[$g]\n";
				}
		}
}
