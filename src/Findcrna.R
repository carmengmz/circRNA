
#MODIFY THIS VARS WITH DESIRED VALUES
####################################################
num_threads = 10
ref_fasta = "hg38.fa"
ref_gtf = "hg38.gtf"
####################################################

bwa index -a bwtsw hg38.fa

bwa mem –T 19 ref.fa reads.fq > aln-se.sam #(for single-end reads)
bwa mem –T 19 ref.fa read1.fq read2.fq > aln-pe.sam #(for paired-end reads)
                                                     
perl CIRI2.pl -T 10 -I aln-pe.sam -O  outfile -F hg38.fa -A hg38.gtf > ciri.log
