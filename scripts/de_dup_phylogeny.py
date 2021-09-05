import pandas as pd
import sys
import numpy as np

df = pd.read_csv(sys.argv[1])
df.replace(float("inf"), np.nan, inplace=True)
df = df.groupby("id").aggregate(max)
df.to_csv("phylogeny.csv")