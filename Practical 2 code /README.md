# This is all the code used in practical code. This focuses only on bash/Linux code.

### Bowtie2:
bowtie2-build ~/practical_2/contigs.fasta newindex 
bowtie2 -x newindex --very-sensitive -1 first_in_pair -2 second_in_pair -S output.sam 

### Convert SAM to BAM:
samtools view -b -o output.bam input.sam 
samtools flagstat output.bam 

### HTSeq:
~/.local/bin/htseq-count -f bam -a 0 -s no -r name -t CDS output.bam prokkagff2gtf.sh > output
