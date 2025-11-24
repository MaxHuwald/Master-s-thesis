#!/bin/bash
#SBATCH --job-name=alignment  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinecpu
#SBATCH --array=1-2
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

NAME=HSV2_333_ARPE19_WDX-3
IN=/project/sysviro/users/Max/analyses/Trimmer_test/
OUT=/project/sysviro/users/Max/analyses/QC
VERSION=1.1.1
GENOME=/project/sysviro/data/reference_genomes/ENO2_RNACS.fasta
BARCODE=/project/sysviro/users/Max/WarpDemuX/$NAME/Demultiplexed
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
### Filter pt > 0 reads and align, keep tags
#samtools view -h $IN/"$NAME".sup-allMods.trimAdapters.dorado.$VERSION.bam | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
#minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
#samtools view -bS - > $OUT/"$NAME".sup-allMods.trimAdapters.dorado."$VERSION".filtered.aligned.bam
### Keep primary aligned only
#samtools view -b -F 2308 $OUT/"$NAME".sup-allMods.trimAdapters.dorado."$VERSION".filtered.aligned.bam > $OUT/"$NAME".sup-allMods.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam



conds=("2h" "CHX")
#barcodes=("barcode04.p0.95.list.txt" "barcode05.p0.95.list.txt")

COND=${conds[$SLURM_ARRAY_TASK_ID-1]}
#BARCODE_LIST=${barcodes[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing condition $COND"

#filterbyname.sh in=$IN/"$NAME".sup-allMods.trimAdapters.dorado.$VERSION.trimmed.bam \
 # out=$BARCODE/"$NAME".sup-allMods."$COND".trimAdapters.dorado.$VERSION.bam \
 # names=$BARCODE/$BARCODE_LIST include=true overwrite=true usejni=t -Xmx16g

samtools view -h $BARCODE/"$NAME".sup-allMods.$COND.trimAdapters.dorado.$VERSION.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
samtools view -b - | \
samtools fastq -Tpt - | \
minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
samtools view -bS - | \
samtools view -b -F2308 - | \
samtools sort -o $OUT/"$NAME".sup-allMods.$COND.trimAdapters.dorado."$VERSION".filtered.aligned.ENO2.primary.bam -

samtools view -b $OUT/"$NAME".sup-allMods.$COND.trimAdapters.dorado."$VERSION".filtered.aligned.ENO2.primary.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/"$NAME".sup-allMods.$COND.trimAdapters.dorado."$VERSION".filtered.aligned.ENO2.primary.bedgraph





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
