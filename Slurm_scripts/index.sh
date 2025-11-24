#!/bin/bash
#SBATCH --job-name=index  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
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

DIREC=/project/sysviro/users/Max/analyses/Conditions/hg38/non_eU

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

# Loop through BAM files and index each
for bamfile in "$DIREC"/*.bam
do
  if [[ -f "$bamfile" ]]; then
    echo "Processing $bamfile ..."
    

    # Index the sorted BAM
    samtools index $bamfile

    echo "Indexed: $bamfile"
  fi
done

echo "All BAM files indexed."

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
