#!/bin/bash

bash code/header2.sh
BASE_DIR=$1
REF_FA=$2
GENOME_SIZE=$3
k=$4
READS_FQ=$5

# Build 'reference'/polished assembly graph and stampy hash
mkdir -p $BASE_DIR/ref/stampy
mkdir -p $BASE_DIR/ref/ctx_bins 
ls $REF_FA > $BASE_DIR/ref.list
echo "cortex ref graph"
cortex_var_31_c1 --kmer_size $k \
--mem_height 20 --mem_width 100 \
--se_list $BASE_DIR/ref.list \
--max_read_len 10000 \
--dump_binary $BASE_DIR/ref/ctx_bins/ref.k$k.ctx \
--sample_id REF
stampy.py -G $BASE_DIR/ref/stampy/REF $REF_FA
stampy.py -g $BASE_DIR/ref/stampy/REF -H $BASE_DIR/ref/stampy/REF

# Make INDEX for read fastqs
# Need to make it work for multiple fastqs again...
rm $BASE_DIR/reads_fq.list
for f in $READS_FQ; do echo -e "$(readlink -f $f)" &>> $BASE_DIR/reads_fq.list; done
echo -e "reads_fq\t$(readlink -f $BASE_DIR/reads_fq.list)\t.\t." &> $BASE_DIR/INDEX

# If running in windows, don't need as much memory
mem_height=16
if [ $GENOME_SIZE -gt 1000000 ]; then
mem_height=23
elif [ $GENOME_SIZE -lt 5001 ]; then
mem_height=14
fi

# Run cortex run_calls with bubble caller
perl ~/apps/cortex/scripts/calling/run_calls.pl --first_kmer $k \
--fastaq_index $BASE_DIR/INDEX \
--manual_override_cleaning $BASE_DIR/CLEANFILE \
--bc yes --pd no \
--outdir $BASE_DIR/results \
--outvcf bubble_caller_k$k --ploidy 1 \
--stampy_hash $BASE_DIR/ref/stampy/REF \
--stampy_bin /apps/well/stampy/1.0.24-py2.7/stampy.py \
--list_ref_fasta $BASE_DIR/ref.list \
--refbindir $BASE_DIR/ref/ctx_bins/ \
--genome_size $GENOME_SIZE --mem_height $mem_height --mem_width 100 \
--vcftools_dir /home/rmnorris/apps/vcftools_0.1.13/ \
--do_union yes --ref CoordinatesAndInCalling --workflow independent \
--logfile $BASE_DIR/log_run_calls_bc_k$k.txt \
--gt_assemblies no

# Update the phillippy + nanopolish assembly
perl code/filter_cortex_vcf.pl $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf
bgzip $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered
tabix -p vcf $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz
cat $REF_FA | vcf-consensus $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_BC_calls_at_all_k.raw.vcf.filtered.gz > $BASE_DIR/ref_bc_k$k.fa 2>$BASE_DIR/log_vcf_consensus_k$k.txt

# Build 'reference'/polished assembly graph and stampy hash of the bc corrected
#ls $BASE_DIR/ref_bc_k$k.fa > $BASE_DIR/ref2.list
#echo "cortex ref graph2"
#cortex_var_31_c1 --kmer_size $k \
#--mem_height 20 --mem_width 100 \
#--se_list $BASE_DIR/ref2.list \
#--max_read_len 10000 \
#--dump_binary $BASE_DIR/ref/ctx_bins/ref2.k$k.ctx \
#--sample_id REF2
#stampy.py -G $BASE_DIR/ref/stampy/REF2 $BASE_DIR/ref_bc_k$k.fa
#stampy.py -g $BASE_DIR/ref/stampy/REF2 -H $BASE_DIR/ref/stampy/REF2

# Run cortex run_calls with path divergence
#echo "pd cortex"
#perl ~/apps/cortex/scripts/calling/run_calls.pl --first_kmer 31 \
#--fastaq_index $BASE_DIR/INDEX \
#--manual_override_cleaning $BASE_DIR/CLEANFILE \
#--bc no --pd yes \
#--outdir $BASE_DIR/results \
#--outvcf path_divergence_k$k --ploidy 1 \
#--stampy_hash $BASE_DIR/ref/stampy/REF2 \
#--stampy_bin /apps/well/stampy/1.0.24-py2.7/stampy.py \
#--list_ref_fasta $BASE_DIR/ref2.list \
#--refbindir $BASE_DIR/ref/ctx_bins/ \
#--genome_size $GENOME_SIZE --mem_height $mem_height --mem_width 100 \
#--vcftools_dir /home/rmnorris/apps/vcftools_0.1.13/ \
#--do_union yes --ref CoordinatesAndInCalling --workflow independent \
#--logfile $BASE_DIR/log_run_calls_pd.txt \
#--gt_assemblies no

# Update the phillippy + nanopolish assembly
#perl code/filter_cortex_vcf.pl $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_PD_calls_at_all_k.raw.vcf
#bgzip $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_PD_calls_at_all_k.raw.vcf.filtered
#tabix -p vcf $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_PD_calls_at_all_k.raw.vcf.filtered.gz
#cat $BASE_DIR/ref_bc_k$k.fa | vcf-consensus $BASE_DIR/results/vcfs/bubble_caller_k"$k"_wk_flow_I_RefCC_FINALcombined_PD_calls_at_all_k.raw.vcf.filtered.gz > $BASE_DIR/ref_pd_k$k.fa 2>$BASE_DIR/log_vcf_consensus_pd_k$k.txt

bash code/footer.sh
