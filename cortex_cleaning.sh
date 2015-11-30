k=15

for sample in JR_FAA56953_22072015_Kox_CAV1374	JR_FAA35548_13042015_KpneCAV1016_clean JR_FAA55103_29072015_CFRE_CAV1741 JR_FAA55652_06082015_kpne_17B15G JR_FAA38229_13042015_KpneCAV1015_clean JR_FAA56052_06082015_ecol_P46212 JR_FAA63668_14102015_kpne_CAV1596 JR_FAA48708_08072015_Ecl_CAV1411 JR_FAA56170_22072015_kpne_CAV1596 JR_FAA48971_08072015_Ecol_P46212 JR_FAA56943_29072015_SMAR_CAV1492; do
echo "/data3/projects/nanopore/gram_negs/initial_eval_2015/nanopore/data/fasta/$sample/all/2d.fasta" > cortex_cleaning/$sample.list

cortex_var_31_c1 --kmer_size $k \
--mem_height 20 --mem_width 100 \
--se_list cortex_cleaning/$sample.list \
--max_read_len 10000 \
--dump_binary cortex_cleaning/$sample.k$k.ctx \
--dump_covg_distribution cortex_cleaning/$sample.k$k.ctx.covg \
--sample_id $sample.q0

cortex_var_31_c1 --kmer_size $k \
--mem_height 20 --mem_width 100 \
--se_list cortex_cleaning/$sample.list \
--max_read_len 10000 \
--quality_score_threshold 10 \
--dump_binary cortex_cleaning/$sample.k$k.q2.ctx \
--dump_covg_distribution cortex_cleaning/$sample.k$k.q2.ctx.covg \
--sample_id $sample.q2

cortex_var_31_c1 --kmer_size $k \
--mem_height 20 --mem_width 100 \
--se_list cortex_cleaning/$sample.list \
--max_read_len 10000 \
--quality_score_threshold 10 \
--dump_binary cortex_cleaning/$sample.k$k.q10.ctx \
--dump_covg_distribution cortex_cleaning/$sample.k$k.q10.ctx.covg \
--sample_id $sample.q10

done &> cortex_cleaning/log_get_covgs.txt
