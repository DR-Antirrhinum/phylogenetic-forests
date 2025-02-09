# phylogenetic-forests
Various scripts used in the genomic analyses carried out in Richardson _et. al_ (2025). 

## Mapping and Pileup of samples
The basic pipeline for preparing data for analysis is as follows:

```mermaid
graph LR;
    FASTQ_1-->BAM_1-->Pileup;
    FASTQ_2-->|BWA| BAM_2-->|Samtools mpileup| Pileup-->|mpileup2sync| SYNC;
    FASTQ_3-->BAM_3-->Pileup;
    
```
### Map reads using BWA-MEM:

    source bwa-0.7.17

    # reference genome
    ref=Am_2019.fasta

    # FASTQ paths
    filepath_f1=$1
    filepath_r2=$2
    in_file_f1=$(basename $filepath_f1)
    in_file_r2=$(basename $filepath_r2)

    outfile=$(echo $in_file_f1 | cut -f 1 -d "." )

    # -M for Picard compatibility
    # -t threads
    # -R Complete read group header line

    EX_READ=$(zcat $in_dir/$in_filepath_f1 | head -n 1)
    ID=$(echo $EX_READ | cut -f3 -d ":")
    FL=$(echo $EX_READ | cut -f4 -d ":")
    RL=$(echo $EX_READ | cut -f10 -d ":")

    bwa mem -M -t 8 -R "@RG\tID:${ID}.LANE${FL}\tSM:${outfile}\tLB:${outfile}\tPL:ILLUMINA\tPU:${ID}.${FL}.${RL}" $ref $in_dir/$filepath_f1 $in_dir/$filepath_r2 > $AlignmentsDir/$outfile.bwa.sam

### Sort SAM files and convert to BAM

    source samtools-1.7

    filepath=$AlignmentsDir/$outfile.bwa.sam

    samtools sort -@ 8 -o $AlignmentsDir/$outfile.bwa.sorted.bam $filepath

### Remove PCR duplicates and carry out local realignment around indels

    source jre-7.21

    ref=Am_2019.fasta 
    pic_path=/picard/1.134/x86_64/jars/
    GATK_path=/GATK/3.5.0/x86_64/jars/

    filepath=$AlignmentsDir/$outfile.bwa.sorted.bam

    #Index sorted bamfile
    samtools index $filepath

    #remove PCR duplicates
    java -Xmx16g -jar ${pic_path}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 TMP_DIR=/tmp INPUT=$filepath OUTPUT=${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.bam  METRICS_FILE=${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.metrics

    samtools index ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.bam

    #local realignment around indels
    java -Xmx16g -jar ${GATK_path}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} -I ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.bam -o ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.realign.intervals

    java -Xmx16g -jar  ${GATK_path}/GenomeAnalysisTK.jar -T IndelRealigner -R ${ref}  -targetIntervals ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.realign.intervals -I ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.bam --out ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.realign.bam

    samtools index ${AlignmentsDir}/${outfile}.bwa.sorted.rmdup.realign.bam

### Create Pileup of BAM files

    source samtools-1.7

    ref=Am_2019.fasta 

    # -q and -Q are read and mapping qualities
    # -t DP to output per-sample read depth
    # -B to disable probabilistic realignment for the computation of base alignment quality (Phred-scaled probability of a read base being misaligned)
    # -A to count orphans (anomalous read pairs in variant calling) 
    # use -r to specify a chromosome / scaffold and submit each separately.

    inChr=$1

    samtools mpileup  -q 40 -Q 30 -t DP -BA -f $ref -r ${inChr} Sample_1.bam Sample_2.bam Sample_n.bam > ${inChr}.pileup


#### Convert Pileup to Popoolation2 SYNC file

See https://github.com/popgenvienna/popoolation2/blob/master/mpileup2sync.jar

    source jre-7.21

    java -Xmx16g -jar mpileup2sync.jar --input ${inChr}.pileup --output ${inChr}.sync --min-qual 30 --threads 1

The output SYNC file is used for treeXY analysis (https://github.com/DR-Antirrhinum/treeXY)

## Generating genetic distance trees from treeXY output files

See _grouping_tree_scan.R_

## Simulating a selective sweep on whole chromosome data

See _artificial_sweep.py_

## Generate whole genome trees with _d<sub>XY</sub>_ and _D_

First, run treeXY_whole_genome.py on the SYNC files of all the chromosomes you wish to analyse. Then, use _WG_dXY_D_trees.R_ to read and generate mean _d<sub>XY</sub>_ and _D_ trees for the whole genome.

## Generate whole genome Maximum Likelihood tree from SYNC files

First, run treeXY with _--write-sync_ enabled to generate treeXY_filtered SYNC files. Next, generate consensus FASTA files for each SYNC file using _consensus_from_sync.py_. Then, run the commands from _process.sh_ to generate a whole genome multi FASTA file for all taxa. This FASTA file can then serve as the input for your phylogenetic software of choice, such as RAxML-NG (Kozlov _et al._, 2019):

    raxml-ng --all --msa WG_treeXY_filtered.fa --model GTR+G --tree pars{10} --bs-trees 100

## Calculating shortest root branch and the cophenetic correlation coefficient for hierarchical clustering trees

See _grouping_tree_scan.R_

## Grouping hierarchical clustering trees based on root division

See _grouping_tree_scan.R_
