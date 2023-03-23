import matplotlib.pyplot as plt
import pandas as pd
import sys 

if len(sys.argv) == 1:
    raise Exception("No argument given")

fname = sys.argv[1]

try:
    df = pd.read_csv(fname)
except Exception as e:
    print(e)
    exit(1)

pairs = df[['n1','n2']].apply(tuple, axis=1)
df = pd.concat([pairs, df[['jac']]], axis=1).rename(columns={0:'pair'})
print(df)

df.plot.bar(x='pair', y='jac')
plt.show()