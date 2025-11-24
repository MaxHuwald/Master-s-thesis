#!/bin/bash
#SBATCH --job-name=modkit_HSV2  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --mem=256gb # Job memory request
#SBATCH --time=06:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinecpu
#SBATCH --cpus-per-task=2
#SBATCH --array=0-5

# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> EXPLANATION
#This script provides mRNA modification workup via modkit
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

#  ------------------------------------------------------------------------------------------------

    echo -e "\n"
    echo "Modules loaded successfully."

    echo -e "\n"
    echo -e "           ============================================"
    echo -e "\n\n"


# /================================================================================================\
# |                                             CODE                                               |
# |________________________________________________________________________________________________|

COND=hg38
EU=eU

# Directory containing BAM files
IN_DIR="/project/sysviro/users/Max/analyses/Conditions/$COND/$EU"

# Output directory
OUT_DIR="/project/sysviro/users/Max/analyses/Modkit/$COND/$EU"
MODKIT="/project/sysviro/bin/modkit"

# Collect all BAM files in the directory
mapfile -d $'\0' files < <(find "$IN_DIR" -maxdepth 1 -type f -name "*.bam" -print0 | sort -z)

# Check that SLURM_ARRAY_TASK_ID is within range
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#files[@]}" ]; then
echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) is out of range (0-${#files[@]})."
exit 1
fi

# Get BAM file for this array job
BAM="${files[$SLURM_ARRAY_TASK_ID]}"
NAME=$(basename "$BAM" .bam)

echo "Processing $BAM with sample name $NAME"

# Run Modkit
"$MODKIT" sample-probs "$BAM" --hist -o "$OUT_DIR" --prefix "${NAME}" \
--only-mapped --percentiles 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9 -f 1.0

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
