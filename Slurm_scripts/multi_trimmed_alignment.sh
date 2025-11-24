#!/bin/bash
#SBATCH --job-name=trim_followup  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinecpu
#SBATCH --array=0-4
# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> EXPLANATION
#Takes Basecalled .bam, filters pt>0, alignes to HSV2 Strain 333 genome, keeps pt, MM and ML flags
#Demultiplexes aligned data set
# (In addition to variables) needs:
#   - FILE 1 [TYPE]
#   - FILE 2 [TYPE]


# /================================================================================================\
# |                                           VARIABLES                                            |
# |________________________________________________________________________________________________|

NAME=Max_HSV2
OUT=/project/sysviro/users/Max/analyses/wiggle_test/
VERSION=1.1.1
IN=/project/sysviro/users/Max/analyses/wiggle_test/"$NAME".sup.trimAdapters.dorado."$VERSION".trimmed_wiggle_${SLURM_ARRAY_TASK_ID}.bam
GENOME=/project/sysviro/data/reference_genomes/HSV2-st333-LS480640.fasta

# /================================================================================================\
# |                                         INITIALISATION                                         |
# |________________________________________________________________________________________________|

    starttime=$(date +"%D %T")
    echo -e "\n"
    echo -e "           ============================================"
    echo -e "           =     started at $starttime         ="
    echo -e "           ============================================"

    echo "<<< Start Job $SLURM_ARRAY_TASK_ID on $HOSTNAME >>>"  # Display job start information
    echo "Script successfully started."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                            MODULES                                             |
# |________________________________________________________________________________________________|

# module load EXAMPLE
module load SAMtools
module load minimap2
module load BBMap 
module load bedtools
#  ------------------------------------------------------------------------------------------------

    echo -e "\n"
    echo "Modules loaded successfully."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                             CODE                                               |
# |________________________________________________________________________________________________|

echo "Processing condition" 

### Filter pt > 0 reads and align, keep tags
#samtools view -h $IN | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
samtools fastq -TXB $IN | \
minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
samtools view -bS - | \
samtools view -b -F2308 - > $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".trimmed_wiggle_${SLURM_ARRAY_TASK_ID}.aligned_HSV2.primary.bam


#samtools view -b -f16 $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam

#samtools view -b -F2324 $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam

#samtools view -b $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bedgraph

#samtools view -b $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bedgraph

#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bed

#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam > $OUT/Visualisation/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bed


# /================================================================================================\
# |                                             ENDING                                             |
# |________________________________________________________________________________________________|

    echo -e "\n"

    endtime=$(date +"%D %T")
    echo -e "\n"
    echo -e "           ============================================"
    echo -e "           =     finished at $endtime        ="
    echo -e "           ============================================"


# =================================================================================================
