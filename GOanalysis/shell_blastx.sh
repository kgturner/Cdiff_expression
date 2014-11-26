################# shell_blast.sh ####################
# Created by Kay Hodgins
# 
#
# Assumes: fasta file
#
# Run: sh shell_blast.sh <DIRECTORY NAME>
#	blastx against the uniprot protein database (tablular output, 1e-10, 20 hits and 10 processors)
#	awk keeps only percent ID > 70 and alignment length >30 sorted by pid
# Reports: blastx results sorted.Argv[1] 
#
########################################################



#! /usr/bin/bash



blastall -p blastx -d  /data/databases/Plant_UniProt -i  $1 -m 8 -e 1e-10 -a 10 -b 20 -o $2

awk '{if ( $4>30) print $1"\t"$2"\t"$3}' $2 > temp 
sort -k 1,1 -k3,3nr temp > sorted.$2



rm temp





