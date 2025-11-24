#!/bin/bash
#SBATCH --job-name=  # Job name
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

BARCODES=/project/sysviro/users/Max/WarpDemuX/HSV2_333_ARPE19_WDX-3/Demultiplexed
CSV=/project/sysviro/users/Max/analyses/HSV2_333_ARPE19_WDX-3/RNAkinet/

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

for txtfile in $BARCODES/barcode04.p0.95.list.txt $BARCODES/barcode05.p0.95.list.txt; do
#$BARCODES/barcode06.p0.95.list.txt $BARCODES/barcode11.p0.95.list.txt; do
    base=$(basename "$txtfile" .txt)
    output="$CSV/${base}.eUfiltered.csv"
    awk '
        NR==FNR { barcodes[$1]; next }           # Read barcodes into array
        FNR==1 { print; next }                    # Print header line from second file
        ($1 in barcodes)                          # Print lines where first field matches barcode
    ' "$txtfile" FS=',' OFS=',' "$CSV/HSV2_333_ARPE19_WDX-3.all.csv" > "$output"
    echo "Filtered output saved to $output"
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
