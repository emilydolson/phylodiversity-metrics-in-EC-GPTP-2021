import glob
import os
import pandas as pd

# This warning is incorrectly being tripped - when dealing
# with time_solved, we take a single-row slice of a data frame
# and then modify that slice. Pandas is worried about those modifications
# not making it back to the original dataframe, but we don't actually care
# if they do.
pd.options.mode.chained_assignment = None  # default='warn'

def process_file(csvfile):
    df = pd.read_csv(csvfile)
    if "generation" not in df:
        df.rename(columns={"update":"generation"}, inplace=True)
    if "generation" not in df:
        return
    # duplicates can happen if problem is solved on time point that's already recorded
    df.drop_duplicates("generation", inplace=True) 

    df.set_index("generation", inplace=True)
    if (csvfile.endswith("systematics.csv")):
        col_dict = {}
        if "phenotype" in csvfile:
            for col in df.columns:
                col_dict[col] = "phenotype_" + col
        else:
            for col in df.columns:
                col_dict[col] = "genotype_" + col

        df.rename(columns=col_dict, inplace=True)

    return df


all_data = []
time_solved_data = []

for path in glob.glob("*/*/[0-9]*"):
    if not os.path.exists(path+"/phylodiversity.csv"):
        continue
    print(path)
    csvs = glob.glob(path+"/"+"*.csv")

    # Doing this as a list comprehension is allegedly faster
    local_dfs = [process_file(csvfile) for csvfile in csvs]

    df = pd.concat(local_dfs, axis=1)    
    local_data = {}

    with open(path + "/run.log") as run_log_file:
        for line in run_log_file:
            if line.startswith("0"):
                break
            elif not line.startswith("set"):
                continue
            line = line.split()
            local_data[line[1]] = line[2]

    for val in local_data:
        df[val] = local_data[val]

    if df.empty:
        continue

    all_data.append(df)

    time_solved = df.index.max()
    is_solved = False
    if os.path.exists(path+"/time_solved"):
        with open(path+"/time_solved") as time_file:
            time_solved = time_file.readlines()[0]
            is_solved = True

    time_solved = int(time_solved)
    time_solved_series = df.loc[time_solved, :]
    time_solved_series["solved"] = is_solved
    time_solved_series["solved_or_finished"] = is_solved or (time_solved == local_data["MAX_GENS"])
    time_solved_data.append(time_solved_series)



res = pd.concat(all_data)
all_time_solved = pd.concat(time_solved_data, axis=1).T
all_time_solved.to_csv("all_time_solved.csv")
res.to_csv("all_data.csv")
