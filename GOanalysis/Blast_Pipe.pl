print "blasting against uniprot proteins and annotating\n";
$file=$ARGV[1];
$dir=$ARGV[0];

#open(IN, "$ARGV[0]"); 


system ("sh shell_blastx.sh $dir$file out.uni_$file");
system ("perl blast_top.pl sorted.out.uni_$file"); #change to shell script to do add file to get top entry
system ("sh shell_join_uni.sh unique.sorted.out.uni_$file");

