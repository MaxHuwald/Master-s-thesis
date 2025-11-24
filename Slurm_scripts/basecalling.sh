#!/bin/bash
#SBATCH --job-name=basecalling.1.1.1  # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=128gb # Job memory request
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinegpu
#SBATCH --gres=gpu:a100-40g:1

# /================================================================================================\
# |                                  REQUIREMENTS & DOCUMENTATION                                  |
# |________________________________________________________________________________________________|

# -> EXPLANATION
# (In addition to variables) needs:
#   - FILE 1 [TYPE]
#   - FILE 2 [TYPE]


# /================================================================================================\
# |                                           VARIABLES                                            |
# |________________________________________________________________________________________________|

IN=/project/sysviro/data/Nanopore-Pod5/Max_HSV2
MODEL=/project/sysviro/users/Max/dorado-models
OUT=/project/sysviro/data/DoradoBasecalling
NAME=HSV2_333_ARPE19_WDX-3
VERSION=1.1.1

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
### Run basecalling
/project/sysviro/dorado-$VERSION-linux-x64/bin/dorado basecaller --trim adapters --estimate-poly-a --models-directory $MODEL/ -r sup,inosine_m6A_2OmeA,pseU_2OmeU $IN/$NAME/ > $OUT/"$NAME".sup-allMods.trimAdapters.dorado.$VERSION.bam

### Generate sequencing summary file 
/project/sysviro/dorado-$VERSION-linux-x64/bin/dorado summary $OUT/"$NAME".sup-allMods.trimAdapters.dorado.$VERSION.bam > $OUT/"$NAME".sup-allMods.trimAdapters.dorado.$VERSION.sequencing_summary.txt
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
