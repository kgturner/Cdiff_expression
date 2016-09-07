################# shell_blast.sh ####################
# Created by Kay Hodgins
# 
#
# Assumes: tabular blast output and /data/databases/ATH_GO_GOSLIM.txt GO terms
#
# Run: sh shell_join.sh <DIRECTORY NAME>
#	
# Reports: blastx results <filename_arabx blastx gets GO terms for the remaining hits <filename>_GO
#	   
#
########################################################

#! /usr/bin/bash



awk '{print $2"\t"$1}' $1 | sed 's/_[A-Z]\+//' > temp


LANG=en_EN sort -k1 temp > tempsort
LANG=en_EN sort -k1 /data/databases/uniprot_go_unique2 > temp2sort
LANG=en_EN join -j1 1 -j2 1 tempsort temp2sort > $1_GO

#rm temp*





