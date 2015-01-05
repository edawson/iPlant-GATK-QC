#!/bin/bash
#SBATCH -J test_gatkqc
#SBATCH -o gatk_test.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4

date

## TACC-specific module commands, if only necessary
module unload samtools
module load python
module load launcher

## Variable for input control in SBATCH go here

numThreads=4
input=bwa_output_WHWT_mem.bam
REFERENCE=umd.tar


# SPLIT_COUNT / 4 = number of records per BWA job
SPLIT_COUNT=4000000
# 8 GB
SORT_RAM=8589934592

tar xzf bin.tgz
PATH=$PATH:`pwd`/bin/
GATK="`pwd`/bin/GenomeAnalysisTK.jar"

qBAM=${input}
qBAMind=${inputIND}
ref=${reference}
output=${output}
performRecal=${RECAL}

#Make temporary folders for splitting
mkdir input
mkdir temp

REFERENCE_F=${REFERENCE} #$(basename ${REFERENCE})
GENOME_F_ARCHIVE=$REFERENCE_F
TARGET="reference"
EXTENSION="amb"
mkdir -p ${TARGET}
case ${GENOME_F_ARCHIVE} in
*.tar.bz2)
    tar xvjf ${GENOME_F_ARCHIVE} -C ${TARGET} ;;
*.zip) 
    unzip -d ${TARGET} -q ${GENOME_F_ARCHIVE} ;;
*.tar.gz)
    tar xvzf ${GENOME_F_ARCHIVE} -C ${TARGET} ;;
*.tgz)
    tar xvzf ${GENOME_F_ARCHIVE} -C ${TARGET} ;;
*.tar)
    echo "The reference is: ${GENOME_F_ARCHIVE}"
    tar xvf ./${GENOME_F_ARCHIVE} -C ${TARGET} ;;
esac

REFERENCE_F=$(basename $(find ${TARGET}/ -name "*.${EXTENSION}" -print0) ".${EXTENSION}")
REFERENCE_F="${TARGET}/$REFERENCE_F"
echo $REFERENCE_F

if [ -n "${alreadyCalled}" ]
    then
    knownCalls="-known ${alreadyCalled}"
else
    knownCalls=""
fi

if [ -n "${LOD}" ]
    then lod="-LOD ${LOD}"
else
    lod=""
fi



if [ -z "${output}" ]
    then output="realigned"
else
    output="${output}.realigned"
fi

## Perform local realignment around indels
## First run RealignmentTargetCreator; this need only be done once
java -Xmx4g -jar ${GATK} \
-T RealignerTargetCreator \
-R ${REFERENCE_F} \
-nt ${numThreads} \
-o ./output.intervals \
${known}

## Cache BAM header for later and split BAM into many sam files.
samtools view -H $qBAM >> bamheader.txt
samtools view $qBAM | split -l $SPLIT_COUNT --numeric-suffixes - input/query

## Add header back to sams and recompress to use with GATK
## Use LaunChair to run in parallel.
rm -rf jobfile.txt
for i in ./input/query*
do
  i=$(basename "${i}")
  echo "cat bamheader.txt ./input/${i} | samtools view -Sb -h - > ./input/reheader_${i}.bam" >> jobfile.txt
done

python ./bin/LaunChair/launcher.py -i jobfile.txt

rm -rf paramlist.aln

count=0
for i in  ./input/reheader_query*
do
  i=$(basename "${i}")
  if [ -n "${performRecal}" ]
  then
    echo "java -Xmx4g -jar ${GATK} -T IndelRealigner -R ${REFERENCE_F} -I ./input/${i} -nt ${numThreads} -targetIntervals ./output.intervals -o ./temp/${output}.${count}.bam ${known} ${lod}; \
          java -Xmx4g -jar ${GATK} -T PrintReads -I ./temp/${output}.${count}.bam -BQSR recal_report.grp -o ./temp/${output}.recal.${count}.bam" >> paramlist.aln
  else
    echo "java -Xmx4g -jar ${GATK} -T IndelRealigner -R ${REFERENCE_F} -I ./input/${i} -nt ${numThreads} -targetIntervals ./output.intervals -o ./temp/${output}.${count}.bam ${known} ${lod}" >> paramlist.aln
  fi
  count=$((count + 1))
done


echo "Launcher...."
date
export TACC_LAUNCHER_SCHED=dynamic
EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE paramlist.aln
date
echo "..Done"


## Merge BAM files back together
echo "Merging sorted BAMs"
OWD=$PWD
cd temp

BAMS=`ls *.sorted.bam`
# Merge and sort
samtools merge ${OWD}/${OUTPUT}.bam ${BAMS} && samtools index ${OWD}/${OUTPUT}.bam
cd $OWD
for i in input temp bin jobfile.txt paramlist.aln
    do
        echo "Removing ${i}"
        rm -rf ${i}
    done
date
