#!/bin/bash
#SBATCH --job-name=post_processing  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=36:00:00 # Time limit hrs:min:sec
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

# -------------------------------
# Input BAM directory (edit or pass as $1)
# -------------------------------
BAM_DIR=/project/sysviro/users/Max/analyses/Conditions/subsampled


# Create output Visualisation folder
VIS_DIR="$BAM_DIR/Visualisation"
mkdir -p "$VIS_DIR"

# -------------------------------
# Loop over BAM files in main folder AND eU/non_eU subfolders
# -------------------------------
find "$BAM_DIR" -maxdepth 2 -type f -name "*.bam" | while read -r bam; do
    base=$(basename "$bam" .bam)   # strip .bam extension

    echo "Processing $base..."

    # Define prefix inside Visualisation folder
    prefix="$VIS_DIR/$base"

    # Reverse strand BAM
    samtools view -b -f16 "$bam" > "${prefix}.rev.bam"

    # Forward strand BAM
    samtools view -b -F2324 "$bam" > "${prefix}.fwd.bam"

    # BedGraph coverage
    samtools view -b "${prefix}.rev.bam" \
      | genomeCoverageBed -ibam stdin -bg -split \
      > "${prefix}.rev.bedgraph"

    samtools view -b "${prefix}.fwd.bam" \
      | genomeCoverageBed -ibam stdin -bg -split \
      > "${prefix}.fwd.bedgraph"

    # BAM → BED12
    bamToBed -bed12 -i "${prefix}.rev.bam" > "${prefix}.rev.bed"
    bamToBed -bed12 -i "${prefix}.fwd.bam" > "${prefix}.fwd.bed"

    echo "  Done: $base"
done

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
