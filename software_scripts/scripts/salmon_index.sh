#!/bin/bash

## assumes star version 2.7.0e

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}
cd ${outpath}

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M29/gencode.vM29.transcripts.fa.gz
zcat gencode.vM29.transcripts.fa.gz |cat - GRCm39.primary_assembly.genome.fa > decoy.aware.gencode.vM29.transcripts.fa
grep "^>" GRCm39.primary_assembly.genome.fa |cut -d " " -f 1 > decoys.txt
sed -i -e 's/>//g' decoys.txt

TP_FASTA="decoy.aware.gencode.vM29.transcripts.fa"
INDEX="salmon_gencode.vM29.index"

module load salmon
call="salmon index -i ${INDEX} -k 31 --gencode -p 8 -t ${TP_FASTA} --decoys decoys.txt"
echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
