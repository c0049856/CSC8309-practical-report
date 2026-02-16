# This document contains all the code from practical 1. Everything was coded in bash/Linux.

### FastQC:
fascqc filename1, filename2

### Assembly:
spades.py -1 filename.gz -2 filename.gz -k 111 --careful --cov-cutoff auto -o contigs

### gnx:
gnx -min 100 -nx 25,50,75 contigs/contigs.fasta

### QUAST:
quast.py contigs.fasta -R FN652779.fa -o quast_results 

### Prokka: 
prokka --outdir ‘STRING’ --kingdom ‘STRING’ --genus ‘STRING’ --mincontiglen 500 
       --usegenus --cpus 2 --locustag ‘PROKKA’ --force contigs.fasta 
