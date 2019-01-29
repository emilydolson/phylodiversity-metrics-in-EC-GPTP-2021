import glob
import os
import pandas as pd

all_data = []

for path in glob.glob("*/*/[0-9]*"):
    if not os.path.exists(path+"/phylodiversity.csv"):
        continue

    csvs = glob.glob(path+"*.csv")

    local_dfs = []
    for csvfile in csvs:
        local_dfs.append(pd.read_csv(csvfile, index_col="generation"))

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

    all_data.append(df)

res = pd.concat(all_data)

res.to_csv("all_data.csv")
