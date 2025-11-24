#!/bin/bash
#SBATCH --job-name=gff3totranscriptome  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=16gb # Job memory request
#SBATCH --time=04:00:00 # Time limit hrs:min:sec
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

DIREC=/project/sysviro/users/Max/Genome
GenePred=/project/sysviro/bin
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

$GenePred/gffread $DIREC/HSV2-st333_Annotation_v1.gff3 -T -o $DIREC/HSV2-st333_Annotation_v1.gtf

#$GenePred/gtfToGenePred $DIREC/HSV2-st333_Annotation_v3.gtf $DIREC/HSV2-st333_Annotation_v3.genePred

#$GenePred/genePredToBed $DIREC/HSV2-st333_Annotation_v3.genePred $DIREC/HSV2-st333_Annotation_v3.bed

#bedtools getfasta -s -split -name \
#    -fi $DIREC/HSV2-st333-LS480640.fasta \
#    -bed $DIREC/HSV2-st333_Annotation_v3.bed \
#    > $DIREC/HSV2-st333_Annotation_v3_transcripts.fa

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
