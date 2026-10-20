"""This script generates SLURM submission scripts for alpha-quartz VASP jobs.

Scripts are stored under:
    /public/home/hkust_jwliu_1/liuxi/script/pyCode/efnn_mc/SiO2/alpha_quartz/
        ├── increase_length/
        ├── decrease_length/
        └── perturb_all/

Each script cd's to the corresponding data directory under:
    /public/home/hkust_jwliu_1/liuxi/Documents/pyCode/efnn_mc/SiO2_data/alpha_quartz/...
"""


from pathlib import Path
# ==============================================================================
# 1. Base Paths & Configuration
# ==============================================================================
# Where the .sh scripts will be written
SCRIPT_ROOT=Path("./scripts_alpha_quartz")
# Where the calculations and VASP input files (INCAR, POSCAR, etc.) reside
DATA_ROOT = Path("/public/home/hkust_jwliu_1/liuxi/Documents/pyCode/efnn_mc/SiO2_data/alpha_quartz")

# Number of configurations per category (matching gen_alpha_quartz_configurations.py)
N_CONFIGS = {
    "increase_length": 5,
    "decrease_length": 5,
    "perturb_all": 5,
}


def gen_slurm_template(job_name, out_log, err_log, work_dir):
    # Must NOT have leading spaces before #SBATCH, otherwise SLURM treats them as regular comments
    return f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p lzicnormal
#SBATCH -N 1                       # Exactly 1 whole node
#SBATCH -n 24                      # Use all 24 cores on the node
#SBATCH --mem=80GB                 # Memory request
#SBATCH -t 0-04:00                 # 4 hours time limit
#SBATCH -o {out_log}
#SBATCH -e {err_log}

# 1. Prevent stack overflow crashes
ulimit -s unlimited

# 2. Load Intel Compiler and Intel MPI
module purge
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0

# 3. Navigate to the working data directory
cd {work_dir}

# 4. Launch calculation
mpirun -np 24 vasp_std
"""


# ==============================================================================
# 2. Top-level Generation (No function wrapper)
# ==============================================================================
for category, n_count in N_CONFIGS.items():
    # Dedicated directory for this category's scripts
    category_script_dir = SCRIPT_ROOT / category
    category_script_dir.mkdir(parents=True, exist_ok=True)
    submission_commands = []
    for ind in range(n_count):
        folder_name = f"{category}_{ind}"
        job_name = f"{category[:4]}_{ind}"  # e.g., incr_0, decr_0, pert_0
        # Target working directory containing POSCAR, INCAR, etc.
        target_data_dir = (DATA_ROOT / folder_name).resolve()
        # Logs are saved directly in the script directory
        out_log = f"out.{folder_name}.out"
        err_log =f"err.{folder_name}.err"
        # Generate script text
        script_content = gen_slurm_template(
            job_name=job_name,
            out_log=str(out_log),
            err_log=str(err_log),
            work_dir=str(target_data_dir),
        )
        # Write individual submission script
        script_file = category_script_dir / f"run_{folder_name}.sh"
        script_file.write_text(script_content)
        submission_commands.append(f"sbatch {script_file.name}")
    # Helper script inside each category folder to submit all its jobs
    submit_helper = category_script_dir / "submit_all.sh"
    submit_helper.write_text("#!/bin/bash\n" + "\n".join(submission_commands) + "\n")
    submit_helper.chmod(0o755)

    print(f"Generated {n_count} scripts in: {category_script_dir}")
