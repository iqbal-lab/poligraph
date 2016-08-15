# poligraph
Polish a nanopore assembly using Illumina reads

1. Take as input the nanopore assembly, and fastq(s) of illumina reads for the assembly.
2. Map the illumina reads to the draft assembly.


The quick way - run poligraph globally:

```perl poligraph.pl --draft_assembly path/to/draft_assembly.fa --reads_fq path/to/reads.fq --base_dir dir_to_work_from --global --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py```

The slower (and better?) way - run in the 10kb windows - right now this is broken  - open issue to fix!

```perl poligraph.pl --draft_assembly path/to/draft_assembly.fa --reads_fq path/to/reads.fq --base_dir dir_to_work_from --cortex_dir path/to/cortex --vcftools_dir path/to/vcftools --stampy_bin path/to/stampy.py```

```
Options:
--draft_assembly path/to/draft_assembly.fa, the assembly to be polished.
--reads_fq path/to/reads.fq, reads to use for polishing. If paired end fastqs, enclose in "".
--base_dir dir_to_work_from, default pwd, this will be treated as the head directory and the poligraph_corrected.fa will be output here. Creates the dir if it doesn't exist.
--window_size int, default 10,000 
--num_procs int, default 20, number of windows to process in parallel. 
--label string, adds a labelled level to directory structure (e.g. so can run with different cleaning methods in the dir and give them a different name).
--cortex_dir path/to/cortex.
--vcftools_dir path/to/vcftools.
--stampy_bin path/to/stampy.py.
--k int, kmer size for de bruijin graph construction, default 31.
--bc yes --pd no, exactly one of bc and pd must be yes. This is the method used by cortex to call variants.
--qthresh int, default 10.```

