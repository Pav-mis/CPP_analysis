import pandas as pd
import sys

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

df = pd.read_table(sys.argv[1], header = None).dropna(axis=0)
groups = df.groupby(by=[3])
for name, group in groups:
    if len(group) > 10: # and name != '-':
        #print(name[0])+'\t'+str(len(group)))
        distances = group[4].values
        components = [[int(x) for x in s.split('.')] for s in distances]
        distanceDf = pd.DataFrame(components)
        protiens = group.iloc[:,0].values
        print(protiens)
        for i in range(1,len(distanceDf.columns)):
            distanceDf[i] = distanceDf[i-1]+distanceDf[i]+1
        length = distanceDf.iloc[:,-1:].values
        distanceDf = distanceDf.div(distanceDf.iloc[:,-1:].values, axis=0)
        distanceDf.insert(0, "length", length)
        distanceDf.insert(0, "names", protiens)
        distanceDf = distanceDf.iloc[:, :-1]
        config = name.replace("-", "_")
        distanceDf.to_csv(config + ".positions.txt", sep="\t", header=False, index=False)