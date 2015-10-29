# poligraph
Polish a nanopore assembly using Illumina reads


1. Map illumina reads to the assembly
2. Take 50kb chunks of the assembly, and pull out the relevant reads
3. For each chunk, call variants against the reference, and fix up the ref.
4. Then take the fixed up ref and run standard WG Cortex calling, and again fix up.
