#ARGV[0] should be full path to directory containing input file ending in '/'
#ARGV[1] should be input file name
#Example usage: perl Arab_Pipe.pl /path/to/input/file/ input.file
#Note: You need to set your database inside shell_blastx_arab.sh as well

print "blasting against a. thaliana proteins and annotating\n";
$file=$ARGV[1];
$dir=$ARGV[0];

#open(IN, "$ARGV[0]"); 
print "Running: sh shell_blastx_arab.sh $dir$file out.ath_$file\n";  
system ("sh shell_blastx_arab.sh $dir$file out.ath_$file"); 
print "Running: perl blast_top.pl sorted.out.ath_$file\n"; 
system ("perl blast_top.pl sorted.out.ath_$file"); #change to shell script to do add file to get top entry
print "Running: sh shell_join.sh unique.sorted.out.ath_$file\n"; 
system ("sh shell_join.sh unique.sorted.out.ath_$file");


