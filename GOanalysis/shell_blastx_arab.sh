################# shell_blast.sh ####################
# Created by Kay Hodgins
# 
#
# Assumes: fasta file
#
# Run: sh shell_blast_arab.sh <DIRECTORY NAME>
#	blastx against the arab protein database (tablular output, 1e-10, 20 hits and 10 processors)
#	awk keeps only  alignment length >30 but can change this depending on what you want to filter by
# Reports: blastx results sorted.ARGV[0] by pid
#
########################################################



#! /usr/bin/bash



blastall -p blastx -d  /home/kgturner/Centaurea_diffusa_expression/GOanalysis/TAIR10_pep_20101214 -i $1 -m 8 -e 1e-10 -a 10 -b 20 -o $2


awk '{if ($4>30) print $1"\t"$2"\t"$3}' $2 > temp 
sort -k 1,1 -k3,3nr temp > sorted.$2



rm temp





