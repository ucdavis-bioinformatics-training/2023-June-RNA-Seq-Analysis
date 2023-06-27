#!/bin/bash

#SBATCH --job-name=salmon_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=60
#SBATCH --mem=15000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=rnaworkshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/salmon-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/salmon-index_%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail@email.com

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}
cd ${outpath}

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz
#zcat gencode.vM32.pc_transcripts.fa.gz |cat - GRCm39.primary_assembly.genome.fa > decoy.aware.gencode.vM32.transcripts.fa
grep "^>" /share/workshop/mrnaseq_workshop/Data/GRCm39.primary_assembly.genome.fa |cut -d " " -f 1 > decoys.txt
sed -i -e 's/>//g' decoys.txt

TP_FASTA="/share/workshop/mrnaseq_workshop/Data/decoy.aware.gencode.vM32.transcripts.fa"
INDEX="salmon_gencode.vM32.index"

module load salmon
call="salmon index -i ${INDEX} -k 31 --gencode -p 8 -t ${TP_FASTA} --decoys decoys.txt"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
