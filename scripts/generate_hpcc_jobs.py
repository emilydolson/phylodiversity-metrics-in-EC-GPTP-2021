import os 
import time
import argparse

def submit_job():
    filename = str(start_seed)+".sb"
    with open(filename, "w") as outfile:
        outfile.write(f"""#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
    
#SBATCH --time=04:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --array={start_seed}-{end_seed}                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=1G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name {name}      # you can give your job a name for easier identification (same as -J)

########## Command Lines to Run ##########
mkdir {dest_dir}
cd {dest_dir}
mkdir $SLURM_ARRAY_TASK_ID
cd $SLURM_ARRAY_TASK_ID

cp {executable} .
./ecology_parameter_sweep -SEED $SLURM_ARRAY_TASK_ID -PROBLEM {problem} -POP_SIZE {pop_size} -SELECTION {selection} -MODES_RESOLUTION 100 -FILTER_LENGTH {pop_size} {additional_args} > run.log
""")
    os.system(f"sbatch {filename}")
    time.sleep(2)

def make_local_job():
    filename = str(start_seed)+".sh"
    with open(filename, "w") as outfile:
        outfile.write(f"""#!/bin/bash
mkdir {dest_dir}
cd {dest_dir}

for seed in {{{start_seed}..{end_seed}}}
do
    echo $seed
    mkdir $seed
    cd $seed

    cp {executable} .
    ./ecology_parameter_sweep -SEED $seed -PROBLEM {problem} -POP_SIZE {pop_size} -SELECTION {selection} -MODES_RESOLUTION 100 -FILTER_LENGTH {pop_size} {additional_args} > run.log
    cd ..
done
""")

# -SELECTION 4 -POP_SIZE 100 -START_POP_SIZE 100 -PROBLEM 4 -N 10 -RESOURCE_SELECT_COST 3 -RESOURCE_SELECT_NICHE_WIDTH .5  -MAX_GENS 5000 -MUT_RATE .01 -RESOURCE_SELECT_FRAC .01 -RESOURCE_SELECT_RES_INFLOW 50 -ECOLOGY_DATA_RES 1

parser = argparse.ArgumentParser()
parser.add_argument('--selections', type=int, nargs='+', help='selection types to use')
parser.add_argument('--problems', type=int, nargs='+', help='problems to use')
parser.add_argument('--additional_args', default="", type=str, help='string containing additional command line args')
parser.add_argument('--pop_size', default=100, type=int, help='population_size')
parser.add_argument('--local', default=False, action="store_true", help='generate local script')
args = parser.parse_args()

additional_args = args.additional_args
start_seed = 0
if os.path.exists(".next_seed"):
    with open(".next_seed") as seedfile:
        start_seed = int(seedfile.readlines()[0])

sel_map = {0:"tournament", 1:"sharing", 2:"lexicase", 3:"ecoea"}
problem_map = {0:"nk", 3:"sorting", 4:"logic", 1:"programsynthesis", 2:"realvalue"}

reps = 10
executable = "/mnt/scratch/dolsonem/ecology_in_EC_parameter_sweep/ecology_parameter_sweep"
if args.local:
    executable = "../../../../ecology_parameter_sweep"
end_seed = start_seed + reps-1
pop_size = args.pop_size

for problem in args.problems:
    for selection in args.selections:
        dest_dir = sel_map[selection] +"/" + problem_map[problem]
        name = sel_map[selection] + problem_map[problem]
        if args.local == True:
            make_local_job()
        else:
            submit_job()
        start_seed = end_seed + 1
        end_seed = start_seed + reps - 1

with open(".next_seed", "w") as seedfile:
    seedfile.write(str(end_seed+1))
