## Running the Python Scripts (Unix / Linux / macOS)

These scripts were custom-developed for the thesis. Follow the steps below to set up a clean Python environment and run them on a Unix-based system.

---

### 1. Create a new virtual environment

Navigate to the directory where you want the environment to live, then create it:

```bash
cd {your_directory}
python3 -m venv {your_environment_name}
```

---

### 2. Activate the environment and install required packages

Activate the virtual environment:

```bash
source {your_environment_name}/bin/activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

*Make sure `requirements.txt` is located in the same directory where you run the command.*

---

### 3. View help information for any script

Move the desired script(s) into your working directory, then display usage instructions:

```bash
python {script_name}.py --help
```

This will show a description of the script and available command-line options.

---

### 4. Run the script

Execute any script using:

```bash
python {script_name}.py [arguments]
```

Provide any required or optional arguments as shown in the `--help` output.


     ### ⚠️ Note on Trimming Script Runtime

The **trimming script** may require **significant runtime** depending on input size and system resources.  
If preferred, you can run it on a compute cluster using SLURM.

An example SLURM script for running the trimmer is available here:

```
Master-s-thesis/SLURM_scripts/Trimmer.sh
```

You may adapt this script to your own cluster environment if needed.
   


   
    
