#!/bin/bash


# Run cgMLST analysis in a loop
#sh cgMLST.sh <inpath>

#Requires parallel

cd $1
for f in *.fna
do
/PATH/TO/lmo-mlst/at_matcher_per_allele.sh $f;

done
