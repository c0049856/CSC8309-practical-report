# This code contains all the additional analysis done to the plasmid in Linux/bash. 
### copying the contigs:
grep -Pzo '(?s)>NODE_3.*>NODE_4' contigs.fasta > contig3.fasta

then nano to remove >NODE_4
### plasmid annotation with prokka: 
prokka --outdir plasmid_annotation --prefix plasmid --genus Chlamydia plasmid.fasta

then blastp to improve annotation -> analysis was moved to RStudio 
