# TagSeq data reduction workflow

## What is TagSeq

3′ TagSeq is a protocol to generate low-cost and exceptionally low-noise gene expression profiling data.

<div id="banner" style="overflow: hidden; display: flex; justify-content:space-around;">
<img src="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/data_reduction/filetypes_figures/tagseq2.png" alt="tagseq_protocol" width="40%"/>
<img src="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/data_reduction/filetypes_figures/tagseq3.png" alt="tagseq_protocol_umi" width="40%" height="90%"/>
</div>

<p align = "right" style="font-family:Times;font-size:12px;">
Weng, etc., Methods Mol Biol. 2022;2398:151-172. doi: 10.1007/978-1-0716-1912-4_13
</p>


The comparison of the traditional RNASeq profile and TagSeq across gene annotations.

<p align = "center">
<img src="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2022-August-RNA-Seq-Analysis/master/data_reduction/filetypes_figures/lexo.png" alt="tagseq_coverage" width="400px"/>
</p>

<p align = "right" style="font-family:Times;font-size:12px;">
https://www.lexogen.com/wp-content/uploads/2015/04/nmeth.f.376.pdf
</p>



### Advantages of TagSeq

  * lower sequencing depth required (~5M vs 20-30M)
  * single-end sequencing is sufficient
  * low cost
  * less sensitivy to sample quality
  * Unique Molecular Identifiers used by default
  * low noise
  * strand specific (>99%)


### Disadvantages of TagSeq

  * require a reference genome with high quality annotation
  * does not have much information on isoform expression
  * only applies to polyadenylated transcripts


### Data reduction workflow

The data reduction steps involved in TagSeq analysis are very similar to regular RNASeq data, except that we should take advantages of the UMIs used to remove any PCR duplicates in the data (not SuperDeduper from htstream). This requires a tool that can identify those UMIs and collapse all duplicates that have the same UMI. [UMItools](https://github.com/CGATOxford/UMI-tools) is one of these tools. It extracts the UMI sequences from raw sequencing reads first. After mapping the reads to the reference genome, it then collapse all identical UMI tagged reads that also map to the same genomic location to one sinlge read. For using with our preprocessing tool htstream seamlessly, we have written our own script to extract the UMI information. After mapping and deduplication, the counts table has to be generated using a separate counting program ([htseq-counts](https://htseq.readthedocs.io/en/release_0.11.1/count.html) or [featurecounts](http://subread.sourceforge.net/)).

1. An example of the preprocessing script used with TagSeq data.

    <pre class="prettyprint"><code class="language-py" style="background-color:333333">#!/bin/bash

    #SBATCH --job-name=salmon_index # Job name
    #SBATCH --nodes=1
    #SBATCH --ntasks=9
    #SBATCH --time=120
    #SBATCH --mem=3000 # Memory pool for all cores (see also --mem-per-cpu)
    #SBATCH --partition=production
    #SBATCH --array=1-38
    #SBATCH --output=slurmout/htstream_%A_%a.out # File to which STDOUT will be written
    #SBATCH --error=slurmout/htstream_%A_%a.err # File to which STDERR will be written

    start=`date +%s`
    echo $HOSTNAME
    echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
    
    sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
    
    inpath="00-RawData"
    outpath="01-HTS_Preproc"
    [[ -d ${outpath} ]] || mkdir ${outpath}
    [[ -d ${outpath}/${sample} ]] || mkdir ${outpath}/${sample}
    
    echo "SAMPLE: ${sample}"
    
    module load htstream/1.3.3
    module load anaconda3/4.9.2
    
    infile=`ls ${inpath} | grep ${sample}_`
    
    echo "FILE: ${infile}"
    
    call="hts_Stats -L ${outpath}/${sample}/${sample}.json -N 'initial stats' \
              -U ${inpath}/${infile} | \
            hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'screen phix' | \
            hts_SeqScreener -A ${outpath}/${sample}/${sample}.json -N 'count the number of rRNA reads'\
            -r -s References/mouse_rrna.fasta | \
            ./extract_UMI_htstream.py --read 1 --length 6 | \
            hts_CutTrim -a 16 -A ${outpath}/${sample}/${sample}.json -N 'trim first 16 bases' | \
            hts_AdapterTrimmer -A ${outpath}/${sample}/${sample}.json -N 'trim adapters' | \
            hts_PolyATTrim --no-left --skip_polyT -A ${outpath}/${sample}/${sample}.json -N 'remove polyA' | \
            hts_QWindowTrim -A ${outpath}/${sample}/${sample}.json -N 'quality trim the ends of reads' | \
            hts_NTrimmer -A ${outpath}/${sample}/${sample}.json -N 'remove any remaining N bases' | \
            hts_LengthFilter -A ${outpath}/${sample}/${sample}.json -N 'remove reads < 50 bp' \
            -n -m 50 | \
            hts_Stats -A ${outpath}/${sample}/${sample}.json -N 'final stats' \
            -f ${outpath}/${sample}/${sample}"

    echo $call
    eval $call
    
    end=`date +%s`
    runtime=$((end-start))
    echo $runtime
    </code></pre>

2. After running this script, you will obtain similar results as we have done in the class, but only for R1 read. Then alignment can be done using STAR (in single-end mode). 

3. Then the deduplication can be performed using UMItools as demonstrated in the script below.


    <pre class="prettyprint"><code class="language-py" style="background-color:333333">#!/bin/bash


    #SBATCH --nodes=1
    #SBATCH --ntasks=8
    #SBATCH --time=360
    #SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
    #SBATCH --partition=production
    #SBATCH --array=1-54
    #SBATCH --output=logs/counts/count_%A_%a.out # File to which STDOUT will be written
    #SBATCH --error=logs/counts/count_%A_%a.err # File to which STDERR will be written
    
    start=`date +%s`
    echo $HOSTNAME
    echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
    
    sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
    annotation="References/gencode.vM26.annotation.gtf"
    inpath='02-STAR_alignment'
    
    echo "SAMPLE: ${sample}"
    
    module load samtools
    
    samtools index -b ${inpath}/${sample}/${sample}_Aligned.sortedByCoord.out.bam
    
    module load umitools
    source activate umitools-1.1.2
    
    umi_tools dedup -I ${inpath}/${sample}/${sample}_Aligned.sortedByCoord.out.bam \
                    --log=${inpath}/${sample}/${sample}_dedup.log \
                    -S ${inpath}/${sample}/${sample}_dedup.bam
    samtools index -b ${inpath}/${sample}/${sample}_dedup.bam
    
    
    module load subread
    
    featureCounts -T 4 --verbose -s 1 \
                  -a ${annotation} \
                  -t exon \
                  -o ${inpath}/${sample}/${sample}_featurecounts.txt \
                  ${inpath}/${sample}/${sample}_dedup.bam \
                  > ${inpath}/${sample}/${sample}_featurecounts.stdout \
                  2> ${inpath}/${sample}/${sample}_featurecounts.stderr
    
    
    end=`date +%s`
    runtime=$((end-start))
    echo $runtime
    </code></pre>

