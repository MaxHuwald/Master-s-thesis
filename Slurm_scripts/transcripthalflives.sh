#!/bin/bash
#SBATCH --job-name=nanopore_RNAkinet # Job name
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=maximilian.huwald@stud.mh-hannover.de # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb # Job memory request
#SBATCH --time=96:00:00 # Time limit hrs:min:sec
#SBATCH --output=/project/sysviro/users/Max/logs/%x_%j.log # Standard output and error log
#SBATCH -p leinegpu
#SBATCH --gres=gpu:a100-40g:1


### Set up variables

# Paths
BAM_DIR="/project/sysviro/users/Max/analyses/Conditions/transcriptome"
CSV_DIR="/project/sysviro/users/Max/analyses/RNAkinet"
OUT_DIR="/project/sysviro/users/Max/analyses/Conditions/halflives"

# Create output dir if needed
mkdir -p "$OUT_DIR"

module load Python/3.10.4-GCCcore-11.3.0
source /project/sysviro/users/Max/python_env/RNAkinet/bin/activate

echo "modules loaded succesfully, environment activated"

# Loop through BAM files
for bam in "$BAM_DIR"/*.bam; do
    # Extract condition tag (XXX) from BAM filename — everything before first '.'
    condition=$(basename "$bam" | cut -d'.' -f1)

    # Skip CHX condition
    if [[ "$condition" == "CHX" ]]; then
        echo "Skipping CHX condition"
        continue
    fi

    # Find matching CSV in predictions directory
    csv="$CSV_DIR/${condition}.p0.95.list.eUfiltered.csv"
    if [[ ! -f "$csv" ]]; then
        echo "No CSV found for condition $condition, skipping"
        continue
    fi


    # Output file name
    out_csv="$OUT_DIR/${condition}.halflives.csv"

    echo "Processing condition $condition (TL=2)"
    rnakinet-predict-halflives \
        --transcriptome-bam "$bam" \
        --predictions "$csv" \
        --tl 2 \
        --output "$out_csv"
done


echo "script finished"

# -------------------------------------------
