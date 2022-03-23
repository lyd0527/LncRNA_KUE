#### Alignment

id=test
db_star=Nitab_star
PE1=test_1.pe
PE2=test_1.pe
nt=30

STAR --twopassMode Basic --outSAMstrandField intronMotif \
               --outFilterIntronMotifs RemoveNoncanonical --outSAMattrIHstart 0 \
               --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --outSAMmapqUnique 60 \
               --genomeDir $db_star --runThreadN $nt \
               --readFilesIn $PE1 $PE2


#############################################################
#### RABT
#############################################################

samtools merge $id1.bam $id2.bam $id3.bam ....> merged.bam
stringtie merged.bam -G nitab.gtf -o Nitab-K.gtf -p 50

#############################################################
#### Quantification
#############################################################

featureCounts -T $nt -p -g gene_id -a Nitab-K.gtf -o Counts.txt *.bam
