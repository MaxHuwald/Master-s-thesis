#!/bin/bash
#SBATCH --job-name=combining  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinecpu

# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> EXPLANATION
# This script provides output files containing transcription start sites (TSS, the 5' end of RNAs) or cleavage & polyadenylation sites (CPAS, the 3'end of RNAs).
# (In addition to variables) needs:
#   - FILE 1 [TYPE]
#   - FILE 2 [TYPE]


# /================================================================================================\
# |                                           VARIABLES                                            |
# |________________________________________________________________________________________________|
trim="trimmed"
DIR="/project/sysviro/users/Max/analyses/Conditions/$trim"

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

#  ------------------------------------------------------------------------------------------------

    echo -e "\n"
    echo "Modules loaded successfully."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                             CODE                                               |
# |________________________________________________________________________________________________|

cd "$DIR" || { echo "Directory not found"; exit 1; }

mkdir -p Visualisation

echo "Merging BAM files in: $DIR"
echo "Excluding: 4h.trimAdapters.dorado.1.1.1."$trim".filtered.aligned.primary.bam"

# Create a list of BAMs to merge (all except the excluded one)
BAMS=$(ls *.bam | grep -v '4h.trimAdapters.dorado.1.1.1.${trim}.filtered.aligned.primary.bam')

echo "-----------------------------------------------------------"
echo "The following BAM files will be merged:"
for f in $BAMS; do
    echo "  - $f"
done
echo "-----------------------------------------------------------"

# Output merged BAM (in Visualisation/)
OUT="${trim}.combined_6samples.bam"

echo "Merging into Visualisation/$OUT"
samtools merge -@ "$SLURM_CPUS_PER_TASK" Visualisation/"$OUT" $BAMS

# Index the merged BAM
samtools index Visualisation/"$OUT"

echo "Finished merging."
echo "Output BAM: Visualisation/$OUT"

########################################
# Subsample merged BAM to 30,000 reads
########################################

TARGET_READS=30000
MERGED_BAM="Visualisation/$OUT"
SUB_BAM="Visualisation/${OUT%.bam}.subsampled_30k.bam"

echo "Subsampling $MERGED_BAM to ~$TARGET_READS reads with seed 42"

# Count total reads in merged BAM
TOTAL_READS=$(samtools view -c "$MERGED_BAM")
echo "Total reads in merged BAM: $TOTAL_READS"

if [ "$TOTAL_READS" -le "$TARGET_READS" ]; then
    echo "Total reads ($TOTAL_READS) <= $TARGET_READS; copying without subsampling."
    cp "$MERGED_BAM" "$SUB_BAM"
    samtools index "$SUB_BAM"
else
    # Compute fraction = 30000 / TOTAL_READS (as floating point)
    FRAC=$(awk -v t="$TARGET_READS" -v n="$TOTAL_READS" \
        'BEGIN{ printf "%.10f", t/n }')

    # samtools -s expects SEED.FRACTION with no leading 0 in the fraction part,
    # so convert 0.00345 -> .00345
    FRAC_NOZERO=$(echo "$FRAC" | sed 's/^0//')

    echo "Sampling fraction: $FRAC (samtools -s 42$FRAC_NOZERO)"

    samtools view -@ "$SLURM_CPUS_PER_TASK" -s 42"$FRAC_NOZERO" -b "$MERGED_BAM" > "$SUB_BAM"

    # Index subsampled BAM
    samtools index "$SUB_BAM"
fi

echo "Subsampled BAM: $SUB_BAM"
echo "Done."

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
