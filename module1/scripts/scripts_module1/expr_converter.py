import GEOparse
import pandas as pd
import numpy as np
import os

probes_conv = GEOparse.parse_GSM(
    "/home/cyril/Documents/Master/sem_1/Case_study/module1/data/GPL6457_old_annotations.txt.gz")

#gse = GEOparse.get_GEO("GSE24616", destdir="./")
gse = GEOparse.get_GEO(filepath="./GSM606890.TXT.GZ.soft")

char = {"stage": [], "time": [], "sex": [], "sample_name": []}
for gsm_name, gsm in sorted(gse.gsms.iteritems()):
    char["stage"].append(gsm.metadata['characteristics_ch1'][1].split(": ")[1])
    char["time"].append(gsm.metadata['characteristics_ch1'][2].split(": ")[1])
    char["sex"].append(gsm.metadata['characteristics_ch1'][3].split(": ")[1])
    char["sample_name"].append(gsm.metadata['title'][0].split(": "))


print(char["stage"][3], char["time"][3],char["sex"][3], char["sample_name"][3])

GPL = gse.gpls["GPL6457"]
pivoted_samples = gse.pivot_samples('VALUE')


strata = pd.read_csv("phylostrata.txt",sep="\t",header=None)
strata.columns = ["GeneID","ProbeID","age"]
strata.set_index("ProbeID",inplace=True)

matched_data = pivoted_samples.join(strata, how="inner").groupby(level=0).last()

unique_data = matched_data.groupby("GeneID").mean()

char_pd = pd.DataFrame(char,index=char["sample_name"])
mixed = char_pd[char_pd.sex == "mixed"].sample_name.tolist()
mixed += char_pd[char_pd.sex == "female"].sample_name.tolist()

char_pd["timing_number"] = 0
time_stamps = char_pd.time.unique()
for i in xrange(len(time_stamps)):
    char_pd.loc[char_pd.time == time_stamps[i],"timing_number"] = i+1

experiment_index = char_pd[char_pd.index.isin(mixed)].reset_index().groupby("timing_number")["index"].apply(lambda x: np.array(x))
set_mean = {}
stages = []
for d, col_list in experiment_index.iteritems():
    set_mean[d] = unique_data[col_list].mean(axis=1)
    stages.append(char_pd[char_pd.index.isin(col_list)].stages[0])

mean_data = pd.DataFrame(set_mean)