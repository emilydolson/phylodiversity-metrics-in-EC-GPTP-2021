import pandas as pd
import sys
import numpy as np

df = pd.read_csv(sys.argv[1])
print(df.groupby("phenotype").count())