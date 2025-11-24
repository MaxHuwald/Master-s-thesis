#!/bin/bash
#SBATCH --job-name=custom_Trimmer  # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%A_%a.log
#SBATCH -p leinecpu
#SBATCH --array=0-4  # Array of jobs for wiggle 0 to 4

# Variables
IN=/project/sysviro/data/DoradoBasecalling/Max_HSV2.sup.trimAdapters.dorado.0.8.0.bam

# Set output file with max_wiggle value from array task ID
OUT=/project/sysviro/users/Max/analyses/wiggle_test/Max_HSV2.sup.trimAdapters.dorado.0.8.0.trimmed_wiggle_${SLURM_ARRAY_TASK_ID}.bam

starttime=$(date +"%D %T")
echo -e "\n"
echo -e "           ============================================"
echo -e "           =     started at $starttime         ="
echo -e "           ============================================"

echo "<<< Start Job $SLURM_ARRAY_TASK_ID on $HOSTNAME >>>"
echo "Script successfully started."
echo -e "\n"
echo -e "           ============================================"
echo -e "\n\n"

module load Python/3.12.3-GCCcore-13.3.0

echo -e "\n"
echo "Modules loaded successfully."
echo -e "\n"
echo -e "           ============================================"
echo -e "\n\n"

# Run trimming with max_wiggle from 0 to 4 (SLURM_ARRAY_TASK_ID variable)
/project/sysviro/users/Max/python_env/Trimmer/bin/python \
   /project/sysviro/users/Max/python_env/Trimmer/Trimmer.py \
        --input $IN \
        --output $OUT \
        --max_wiggle $SLURM_ARRAY_TASK_ID

echo -e "\n"

endtime=$(date +"%D %T")
echo -e "\n"
echo -e "           ============================================"
echo -e "           =     finished at $endtime        ="
echo -e "           ============================================"

