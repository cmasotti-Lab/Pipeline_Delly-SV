#!/usr/bin/bash

#Pipeline atualizado em 25-08-2023
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado

# PARAMETROS OBRIGATORIOS
SCRATCH60="/home/scratch60/vlira_21set2023/"

# DATA=$(date "+%F") # EDITE SE QUISER USAR UMA PASTA DE UMA DATA ESPECIFICA 
DATA='2023-10-28'
OUTPUT_DIR=$SCRATCH60"/Result_Delly_SV."$DATA

INPUT_DIR="/home/scratch60/rtorreglosa_12novl2023/preprocessing_READ_result/"
BAM_FILES=$(find "$INPUT_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam')
JOBS=5

#TOOLS e DATABASES
REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/"
ANNOVAR="$SCRATCH60/tools/annovar/table_annovar.pl"
ANNOVAR_DB="$SCRATCH60/humandb/"

mkdir $OUTPUT_DIR

find "$INPUT_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam' | grep -Pv "ROP-25-|ROP-26-|ROP-27-|ROP-29-" > $OUTPUT_DIR/samples.list
head -3 $OUTPUT_DIR/samples.list >  $OUTPUT_DIR/TOY.samples.list
OUTPUT_LOG="$OUTPUT_DIR.log"

export OUTPUT_DIR
export OUTPUT_LOG
export REF_FASTA
export ANNOVAR
export ANNOVAR_DB

step1_DellyCall (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step1_DellyCall para Amostra: "$NAME" <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  delly call \
    -t=
    -o $OUTPUT_DIR/step1_DellyCall/$NAME.bcf \
    -g $REF_FASTA/Homo_sapiens_assembly38.fasta \
    $SAMPLE 2> $OUTPUT_DIR/step1_DellyCall/$NAME.log
}
export -f step1_DellyCall


step2_DellyMerge (){
  local SAMPLE=$1 
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step2_DellyMerge <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG
  local SAMPLE_BCF=$(find "$OUTPUT_DIR"/step1_DellyCall/ -maxdepth 1 -mindepth 1  -name '*.bcf')

  delly merge \
    -o $OUTPUT_DIR/step2_DellyMerge/ROP-sites.bcf \
    ${SAMPLE_BCF}  2> $OUTPUT_DIR/step2_DellyMerge/ROP-sites.log
}
export -f step2_DellyMerge


step3_DellyCallGenotype (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step3_DellyCallGenotype para Amostra: "$NAME" <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  delly call \
      -o $OUTPUT_DIR/step3_DellyCallGenotype/$NAME.geno.bcf \
      -v $OUTPUT_DIR/step2_DellyMerge/ROP-sites.bcf \
      -g $REF_FASTA/Homo_sapiens_assembly38.fasta \
      $SAMPLE 2> $OUTPUT_DIR/step3_DellyCallGenotype/$NAME.log
}
export -f step3_DellyCallGenotype


step4_bcftoolsMergeGenotype (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step4_bcftoolsMergeGenotype <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG
  local SAMPLE_BCF=$(find "$OUTPUT_DIR"/step3_DellyCallGenotype/ -maxdepth 1 -mindepth 1  -name '*.geno.bcf')

  bcftools merge \
    -m id \
    -O b \
    -o $OUTPUT_DIR/step4_bcftoolsMergeGenotype/ROP-MergeGeno.bcf \
    ${SAMPLE_BCF} 2> $OUTPUT_DIR/step4_bcftoolsMergeGenotype/ROP-MergeGeno.log

  bcftools index $OUTPUT_DIR/step4_bcftoolsMergeGenotype/ROP-MergeGeno.bcf 2>> $OUTPUT_DIR/step4_bcftoolsMergeGenotype/ROP-MergeGeno.log
}
export -f step4_bcftoolsMergeGenotype

step5_DellyFilter (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step5_DellyFilter <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  delly filter \
    -f germline \
    -p \
    --altaf=0.1 \
    --gq=15 \
    -o $OUTPUT_DIR/step5_DellyFilter/ROP-Filter.bcf \
    $OUTPUT_DIR/step4_bcftoolsMergeGenotype/ROP-MergeGeno.bcf 2> $OUTPUT_DIR/step5_DellyFilter/ROP-Filter.log

  bcftools view $OUTPUT_DIR/step5_DellyFilter/ROP-Filter.bcf > $OUTPUT_DIR/step5_DellyFilter/ROP-Filter.vcf
}
export -f step5_DellyFilter


step6_annovar (){
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step6_annovar  <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

    $ANNOVAR  \
      --vcfinput $OUTPUT_DIR/step5_DellyFilter/ROP-Filter.vcf \
      $ANNOVAR_DB --buildver hg38 --remove \
      --protocol refGene,dgvMerged  \
      --operation g,r --arg '-splicing 5', --polish \
      --otherinfo --thread ${JOBS} --outfile $OUTPUT_DIR/step6_annovar/ROP-annovar > $OUTPUT_DIR/step6_annovar/ROP-annovar.log 2> $OUTPUT_DIR/step6_annovar/ROP-annovar.log2

   sed 's/\\x3b/;/g' $OUTPUT_DIR/step6_annovar/ROP-annovar.hg38_multianno.vcf| sed 's/\\x3d/=/g' > $OUTPUT_DIR/step6_annovar/ROP-annovar.hg38_multianno.correct.vcf 

   bcftools query -l $OUTPUT_DIR/step6_annovar/ROP-annovar.hg38_multianno.correct.vcf  > $OUTPUT_DIR/step6_annovar/sampleID.list
}
export -f step6_annovar



echo "                           >>>>>> Starting Pipeline to Run Pipeline_Delly_SV.sh <<<<<<" >> $OUTPUT_LOG
date >> $OUTPUT_LOG

mkdir $OUTPUT_DIR/step1_DellyCall/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step1_DellyCall  "$@"' 'step1_DellyCall'

mkdir $OUTPUT_DIR/step2_DellyMerge/
#step2_DellyMerge

mkdir $OUTPUT_DIR/step3_DellyCallGenotype/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step3_DellyCallGenotype  "$@"' 'step3_DellyCallGenotype'

mkdir $OUTPUT_DIR/step4_bcftoolsMergeGenotype/
#step4_bcftoolsMergeGenotype

mkdir $OUTPUT_DIR/step5_DellyFilter/
step5_DellyFilter

mkdir $OUTPUT_DIR/step6_annovar/
step6_annovar


echo "" >> $OUTPUT_LOG
echo "                           >>>>>> End Pipeline <<< " >> $OUTPUT_LOG
date >> $OUTPUT_LOG
echo "" >> $OUTPUT_LOG
