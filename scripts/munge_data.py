import glob
import os
import pandas as pd

all_data = []
time_solved_data = []

for path in glob.glob("*/*/[0-9]*"):
    if not os.path.exists(path+"/phylodiversity.csv"):
        continue
    print(path)
    csvs = glob.glob(path+"/"+"*.csv")

    local_dfs = []

    for csvfile in csvs:
        local_dfs.append(pd.read_csv(csvfile))
        if "generation" not in local_dfs[-1]:
            local_dfs[-1].rename(columns={"update":"generation"}, inplace=True)
        if "generation" not in local_dfs[-1]:
            local_dfs.pop()
            continue
        local_dfs[-1].set_index("generation", inplace=True)
        if (csvfile.endswith("species_ecology.csv")):
            col_dict = {}
            for col in local_dfs[-1].columns:
                col_dict[col] = "species_" + col
            local_dfs[-1].rename(columns=col_dict, inplace=True)

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

    time_solved = max(df.index)
    is_solved = False
    if os.path.exists(path+"/time_solved"):
        with open(path+"/time_solved") as time_file:
            time_solved = time_file.readlines()[0]
            is_solved = True

    time_solved = int(time_solved)
    time_solved_series = df.loc[time_solved, :]
    time_solved_series["solved"] = is_solved
    time_solved_series["solved_or_finished"] = is_solved or (time_solved == local_data["MAX_GENS"])
    time_solved_data.append()

    all_data.append(df)

res = pd.concat(all_data)
all_time_solved = pd.concat(time_solved_data, axis=1).T
all_time_solved.to_csv("all_time_solved.csv")
res.to_csv("all_data.csv")
