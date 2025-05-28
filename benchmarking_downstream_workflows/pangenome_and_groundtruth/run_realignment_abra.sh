

threads=48
GenomeAnalysisTK="GenomeAnalysisTK.jar"
abra=abra2-2.19.jar
picard="picard2.27.5.jar"


sample="HG001.BBBC20_HGSVC3.L50.chm13.L50.g70"

ref38="ref.fa"
ref38_dic="ref.dict"
samtools faidx ${ref38}
samtools dict ${ref38} > ${ref38_dic}


samtools view -@${threads} -hb -F 4  ${sample}.bam > ${sample}_mapped.bam
samtools index -@${threads}  ${sample}_mapped.bam

cat ${sample}_mapped.bam | bamleftalign --fasta-reference ${ref38} --compressed  > ${sample}.left_shifted.bam 

#samtools addreplacerg -@${threads}  -r "@RG\tID:grch38\tSM:HG002\tPL:ILLUMINA\tDS:novaseq\tPU:novaseq" -w   ${sample}.left_shifted.bam  -o ${sample}.left_shifted_RG.bam 
#samtools index  -@${threads}   ${sample}.left_shifted_RG.bam


java -jar ${picard} ReorderSam -I ${sample}.left_shifted.bam -O ${sample}.left_shifted_reor.bam   --SEQUENCE_DICTIONARY ${ref38_dic} -CREATE_INDEX TRUE -VALIDATION_STRINGENCY LENIENT
java -jar ${GenomeAnalysisTK} -T RealignerTargetCreator --remove_program_records  -drf DuplicateRead    --disable_bam_indexing   -nt ${threads}  -R ${ref38} -I ${sample}.left_shifted_reor.bam  --out forIndelRealigner.intervals 
awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' forIndelRealigner.intervals > forIndelRealigner.intervals.bed

in_expansion_bases=160
bedtools slop -i forIndelRealigner.intervals.bed -g ${ref38}".fai" -b "${in_expansion_bases}" > out.intervals.bed

java -Xmx140G -jar $abra  --targets out.intervals.bed  --in ${sample}.left_shifted_reor.bam   --out ${sample}_abra.bam --ref ${ref38} --index --threads ${threads} 

samtools index  -@ ${threads}  ${sample}_abra.bam 



