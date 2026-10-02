import os
import subprocess
from datetime import datetime, timedelta

def get_sofar_hpss(
    start_date_str, 
    end_date_str, 
    slurm_params=None,
    work_dir="/scratch3/NCEPDEV/marine/Edwin.Rainville/sofar/data/bufr_data"
):
    """
    Generates and submits a Slurm batch script to retrieve and extract daily Sofar 
    buoy data from the NCEP HPSS archive, pre-staging files in chunks via HSI.
    """
    if slurm_params is None:
        slurm_params = {
            "ntasks": 1,
            "queue": "batch",
            "partition": "u1-service",
            "walltime": "01:00:00",
            "account": "marine-cpu",
            "job_name": "get_sofar_hpss"
        }

    archive_base = "/NCEPPROD/hpssprod/runhistory"
    os.makedirs(work_dir, exist_ok=True)
    
    start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.strptime(end_date_str, "%Y-%m-%d")
    
    # Gather all archive paths and handle extraction commands
    archive_files = []
    extraction_commands = []
    current = start_date
    while current <= end_date:
        yy, mm, dd = current.strftime("%Y"), current.strftime("%m"), current.strftime("%d")
        file_atm = f"{archive_base}/rh{yy}/{yy}{mm}/{yy}{mm}{dd}/dcom_{yy}{mm}{dd}.tar"
        target_file = f"sofar_b001_xx121_{yy}{mm}{dd}"
        
        archive_files.append(file_atm)
        
        extraction_commands.append(f"""
target_file="{target_file}"
if [ -f "$target_file" ]; then
    echo "Skipping {yy}-{mm}-{dd}: $target_file already exists."
else
    echo "Processing {yy}-{mm}-{dd}..."
    htar -P -xvf {file_atm} ./b001
    if [ -f "b001/xx121" ]; then
        mv b001/xx121 "$target_file"
        echo "Successfully extracted $target_file"
    else
        echo "Warning: Target xx121 not found in {file_atm}"
        missing_files+=("$target_file")
    fi
    rm -rf b001
fi
""")
        current += timedelta(days=1)

    # Dynamically build the SBATCH header using the dictionary
    job_name = slurm_params.get("job_name", "get_sofar_hpss")
    script_content = f"""#!/bin/sh --login
#SBATCH --ntasks={slurm_params.get("ntasks", 1)}
#SBATCH -q {slurm_params.get("queue", "batch")}
#SBATCH --partition={slurm_params.get("partition", "u1-service")}
#SBATCH -t {slurm_params.get("walltime", "01:00:00")}
#SBATCH -A {slurm_params.get("account", "marine-cpu")}
#SBATCH -J {job_name}
#SBATCH -o {work_dir}/{job_name}_%j.out

cd {work_dir}

# Bash array to track files that fail extraction
missing_files=()

{"".join(extraction_commands)}

echo "All requested dates processed."
echo ""
if [ ${{#missing_files[@]}} -ne 0 ]; then
    echo "=========================================="
    echo "The following files could NOT be extracted:"
    for file in "${{missing_files[@]}}"; do
        echo " - $file"
    done
    echo "=========================================="
else
    echo "All requested files were successfully extracted or already present."
fi
"""

    script_path = os.path.join(work_dir, "submit_hpss.sh")
    with open(script_path, "w") as f:
        f.write(script_content)

    # 4. Submit the job to the service node
    result = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
    print(f"Job Submission Output: {result.stdout.strip()}")
    if result.stderr:
        print(f"Errors: {result.stderr.strip()}")