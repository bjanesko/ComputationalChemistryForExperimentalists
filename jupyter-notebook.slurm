#!/bin/bash -l
#SBATCH --job-name=jupyter-notebook
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --time=01:00:00
#SBATCH --partition=research
#SBATCH --output=jupyter-notebook_%j.out
#SBATCH --error=jupyter-notebook_%j.err

# Get SSH tunneling info
XDG_RUNTIME_DIR=""
node=$(hostname -s)
user=$(whoami)

# Print tunneling instructions into output
echo -e "
# Note: below 8888 is used to signify the port.
#       However, it may be another number if 8888 is in use.
#       Check jupyter_notebook_%j.err to find the port.

# Command to create SSH tunnel:
ssh -N -f -L 8888:${node}:8888 ${user}@hpclogin.it.tcu.edu

# Use a browser on your local machine to go to:
http://localhost:8888/
"

# Launch the Jupyter session
jupyter-notebook --no-browser --ip=${node}

# Sleep an hour until the job is done
sleep 3600
