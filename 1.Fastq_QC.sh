### Cleaning
nt=30
sample_1=test_1.fq.gz
sample_2=test_2.fq.gz

fastp -q 20 -u 40 -M 0 -n 5 -l 50 -w $nt -i $sample_1 -I $sample_2 -o $sample_1.pe -O $sample_2.pe

### QC
./FastQC/fastqc -t 50 *.fq.gz *.pe
