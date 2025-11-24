#!/bin/bash
#SBATCH --job-name=nanopore_RNAkinet # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb # Job memory request
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinegpu
#SBATCH --gres=gpu:a100-40g:1


### Set up variables
IN=/project/sysviro/data/Nanopore-Pod5/Xianming/HFF_mock_HUVEC_mock_wdx_p2-1/20251029_1258_P2S-01949-B_PBE05369_2c641947/pod5
NAME=Mock
OUT=/project/sysviro/users/Max/analyses/$NAME/RNAkinet

# ------------------------------
echo "task successfully started"
# ------------------------------

#RNAkinet
module load Python/3.10.4-GCCcore-11.3.0
source /project/sysviro/users/Max/python_env/RNAkinet/bin/activate

echo "modules loaded succesfully, environment activated"

rnakinet-inference --path $IN --kit r10 --output $OUT/"$NAME".all.csv --format pod5

echo "script finished"

# -------------------------------------------
