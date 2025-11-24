#!/bin/bash
#SBATCH --job-name=eU_bam_transcriptome  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=32gb # Job memory request
#SBATCH --time=08:00:00 # Time limit hrs:min:sec
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

COND=transcriptome
BAM_DIR="/project/sysviro/users/Max/analyses/Conditions/$COND"
CSV_DIR="/project/sysviro/users/Max/analyses/RNAkinet"
OUT_DIR="/project/sysviro/users/Max/analyses/Conditions"
PASS=0.9
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

for bam in "$BAM_DIR"/*.bam; do
base=$(basename "$bam")
base=${base%%.trimAdapters*}
csv="$CSV_DIR/$base.p0.95.list.eUfiltered.csv"

echo "Checking BAM: $bam"
echo "  Expected CSV: $csv"


if [[ -f "$csv" ]]; then
echo "Processing $base..."

# extract passing read IDs
awk -F, -v pass="$PASS" 'NR>1 && $2 >= pass {print $1}' "$csv" > "$OUT_DIR/read_ids_per_eU/$COND/$base.$COND.$PASS.passIDs.txt"

# filter BAM: keep reads in passIDs.txt, write others to a separate BAM
samtools view -b -N "$OUT_DIR/read_ids_per_eU/$COND/$base.$COND.$PASS.passIDs.txt" -U "$OUT_DIR/$COND/non_eU/$base.$COND.non_eU.$PASS.filtered.bam" "$bam" \
> "$OUT_DIR/$COND/eU/$base.$COND.eU.$PASS.filtered.bam"

else
  echo "No matching CSV for $base"
fi
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
