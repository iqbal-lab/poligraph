#!/bin/bash
k=31
label='manual_clean'

sample=JR_FAA63658_29092015_ecol_P46212
ILLUMINA_READS="/data3/projects/nanopore/gram_negs/initial_eval_2015/illumina/P46212/someseq/C00006949_R00000042.bam"
READS="/data3/projects/nanopore/gram_negs/initial_eval_2015/nanopore/data/fasta/$sample/all/2d.fasta"
ICLEANFILE="reads_fq\t31\t5\n"
NCLEANFILE="reads_fq\t31\t5\n"

## PBcR_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/assembly.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_nanopolish#################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/polished.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_2nd_nanopolish#################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/polished2.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/$sample/utg.fa

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_and_nanopolish##################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/$sample/polished.fa

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"


sample=JR_FAA63668_14102015_kpne_CAV1596
ILLUMINA_READS="/data3/projects/nanopore/gram_negs/initial_eval_2015/illumina/CAV1596/miseq/CAV-1596_R00000049_v3.bam"
READS="/data3/projects/nanopore/gram_negs/initial_eval_2015/nanopore/data/fasta/$sample/pass/2d.fasta"
ICLEANFILE="reads_fq\t31\t35\n"
NCLEANFILE="reads_fq\t31\t50\n"

## PBcR_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/assembly.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_nanopolish#################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/polished.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_2nd_nanopolish#################
#DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/$sample/polished2.fa

#DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_illumina_correction/
#mkdir -p $DIR
#echo -e $ICLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

#DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_nanopore_all2d_correction/
#mkdir -p $DIR
#echo -e $NCLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/$sample/utg.fa

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_and_nanopolish##################
#DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/$sample/polished.fa

#DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_illumina_correction/
#mkdir -p $DIR
#echo -e $ICLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

#DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_nanopore_all2d_correction/
#mkdir -p $DIR
#echo -e $NCLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

sample=Ecoli_K12
ILLUMINA_READS="\"/data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R1_001_val_1.fq /data3/projects/nanopore/gram_negs/initial_eval_2015/miseq/ecoli_k12_mg1655/MiSeq/1_S1_L001_R2_001_val_2.fq\""
READS="/data2/users/rachel/projects/initial/loman_method/raw.reads.fasta"
ICLEANFILE="reads_fq\t31\t20\n"
NCLEANFILE="reads_fq\t31\t20\n"

## PBcR_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/sampleDataOxford/assembly.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_nanopolish#################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/sampleDataOxford/polished.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_nanopolish_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## PBcR_and_2nd_nanopolish#################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/phillippy_method/sampleDataOxford/polished2.fa

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/PBcR_and_2nd_nanopolish_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_only###########################
DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/sampleDataOxford/utg.fa

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_illumina_correction/
mkdir -p $DIR
echo -e $ICLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_only_nanopore_all2d_correction/
mkdir -p $DIR
echo -e $NCLEANFILE > $DIR/CLEANFILE
echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"

## miniasm_and_nanopolish##################
#DRAFT_ASSEMBLY=/data2/users/rachel/projects/initial/miniasm_method/sampleDataOxford/polished.fa

#DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_illumina_correction/
#mkdir -p $DIR
#echo -e $ICLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $ILLUMINA_READS &>> $DIR/log_cortex_correction_k$k$label.txt"

#DIR=phillippy_plus_cortex_method/$sample/$label/miniasm_and_nanopolished_nanopore_all2d_correction/
#mkdir -p $DIR
#echo -e $NCLEANFILE > $DIR/CLEANFILE
#echo "bash code/cortex_correction.sh $DIR $DRAFT_ASSEMBLY 4610000 $k $READS &>> $DIR/log_cortex_correction_k$k$label.txt"
