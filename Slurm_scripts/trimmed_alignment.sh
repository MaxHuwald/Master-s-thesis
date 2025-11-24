#!/bin/bash
#SBATCH --job-name=trim_followup  # Job name
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

IN=/project/sysviro/users/Max/analyses/Trimmer_test/Max_HSV2.sup.trimAdapters.dorado.1.1.1.trimmed.bam
NAME=HSV2_333_ARPE19_WDX-2
OUT=/project/sysviro/users/Max/analyses/Trimmer_test/
VERSION=1.1.1
GENOME=/project/sysviro/data/reference_genomes/HSV2-st333-LS480640.fasta
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

echo "Processing condition" 

conds=("6h" "8h")
barcodes=("barcode04.p0.9.list.txt" "barcode05.p0.9.list.txt")

COND=${conds[$SLURM_ARRAY_TASK_ID-1]}
BARCODE_LIST=${barcodes[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing condition $COND"

filterbyname.sh in=$IN \
  out=$BARCODE/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.bam \
  names=$BARCODE/$BARCODE_LIST include=true overwrite=true usejni=t -Xmx16g


### Filter pt > 0 reads and align, keep tags
samtools view -h $BARCODE/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
samtools view -b - | \
samtools fastq -Tpt,XB - | \
minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
samtools view -bS - | \
samtools view -b -F2308 - > $OUT/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.bam


samtools view -b -f16 $OUT/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.rev.bam

samtools view -b -F2324 $OUT/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.fwd.bam

samtools view -b $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.rev.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.rev.bedgraph

samtools view -b $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.fwd.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.fwd.bedgraph

bamToBed -bed12 -i $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.rev.bam > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.rev.bed

bamToBed -bed12 -i $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.fwd.bam > $OUT/Visualisation/"$NAME"."$COND".sup.trimAdapters.dorado."$VERSION".trimmed.filtered.aligned.primary.fwd.bed




### Demultiplex data
#filterbyname.sh in=$IN/"$NAME".sup.trimAdapters.dorado.$VERSION.bam out=$BARCODE/"$NAME".sup.eU.trimAdapters.dorado.$VERSION.bam names=$BARCODE/barcode05.p0.9.list.txt include=true overwrite=true usejni=t -Xmx16g
### Filter pt > 0 reads and align, keep tags
#samtools view -h $BARCODE/"$NAME".sup.eU.trimAdapters.dorado.$VERSION.bam | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
#minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
#samtools view -bS - > $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.bam
### Keep primary aligned only
#samtools view -b -F 2308 $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.bam > $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam
### Create fwd and rev
# Get reverse strand reads only.
#samtools view -b -f16 $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam
# Get forward strand reads only
#samtools view -b -F2324 $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam
# Generate bedgraph files
#samtools view -b $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bedgraph
#samtools view -b $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bedgraph

# Generate bed files
#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bed

#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam > $OUT/Visualisation/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bed




#filterbyname.sh in=$IN/"$NAME".sup.trimAdapters.dorado.$VERSION.bam out=$BARCODE/"$NAME".sup.DMSO.trimAdapters.dorado.$VERSION.bam names=$BARCODE/barcode06.p0.9.list.txt include=true overwrite=true usejni=t -Xmx16g
#samtools view -h $BARCODE/"$NAME".sup.DMSO.trimAdapters.dorado.$VERSION.bam | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
#minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
#samtools view -bS - > $OUT/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.bam
### Keep primary aligned only
#samtools view -b -F 2308 $OUT/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.bam > $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam
### Create fwd and rev
# Get reverse strand reads only.
#samtools view -b -f16 $OUT/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam
# Get forward strand reads only
#samtools view -b -F2324 $OUT/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam
# Generate bedgraph files
#samtools view -b $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bedgraph
#samtools view -b $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bedgraph

# Generate bed files
#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bed

#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam > $OUT/Visualisation/"$NAME".sup.DMSO.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bed




#filterbyname.sh in=$IN/"$NAME".sup.trimAdapters.dorado.$VERSION.bam out=$BARCODE/"$NAME".sup.STM2457.trimAdapters.dorado.$VERSION.bam names=$BARCODE/barcode11.p0.9.list.txt include=true overwrite=true usejni=t -Xmx16g
#samtools view -h $BARCODE/"$NAME".sup.STM2457.trimAdapters.dorado.$VERSION.bam | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
#minimap2 -ax splice -k14 -uf -y --secondary=no $GENOME - | \
#samtools view -bS - > $OUT/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.bam
### Keep primary aligned only
#samtools view -b -F 2308 $OUT/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.bam > $OUT/"$NAME".sup.eU.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam
### Create fwd and rev
# Get reverse strand reads only.
#samtools view -b -f16 $OUT/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam
# Get forward strand reads only
#samtools view -b -F2324 $OUT/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.bam > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam
# Generate bedgraph files
#samtools view -b $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bedgraph
#samtools view -b $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam | genomeCoverageBed -ibam stdin -bg -split > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bedgraph

# Generate bed files
#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bam > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.rev.bed

#bamToBed -bed12 -i $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bam > $OUT/Visualisation/"$NAME".sup.STM2457.trimAdapters.dorado."$VERSION".filtered.aligned.primary.fwd.bed



###Human mapping
#samtools view -h $IN/"$NAME".sup.trimAdapters.dorado.$VERSION.bam | \
#awk 'BEGIN{FS="\t";OFS="\t"} /^@/ {print; next} {pt=-1; for(i=12;i<=NF;i++){if($i ~ /^pt:i:/){split($i,a,":");pt=a[3]}} if(pt>0) print}' | \
#samtools view -b - | \
#samtools fastq -Tpt,MM,ML - | \
#minimap2 -ax splice -k14 -uf -y --secondary=no /project/sysviro/data/reference_genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa - | \
#samtools view -bS - > $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.hg38.bam
### Keep primary aligned only
#samtools view -b -F 2308 $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.hg38.bam > $OUT/"$NAME".sup.trimAdapters.dorado."$VERSION".filtered.aligned.hg38.primary.bam

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
