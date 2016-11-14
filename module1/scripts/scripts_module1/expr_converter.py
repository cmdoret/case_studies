import GEOparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


probes_conv = GEOparse.parse_GSM(
    "/home/cyril/Documents/Master/sem_1/Case_study/module1/data/GPL6457_old_annotations.txt.gz")


gse = GEOparse.get_GEO("GSE24616", destdir="./")
# gse = GEOparse.get_GEO(filepath="./GSM606890.TXT.GZ.soft")

char = {"stage": [], "time": [], "sex": [], "sample_name": []}
for gsm_name, gsm in sorted(gse.gsms.iteritems()):
    char["stage"].append(gsm.metadata['characteristics_ch1'][1].split(": ")[1])
    char["time"].append(gsm.metadata['characteristics_ch1'][2].split(": ")[1])
    char["sex"].append(gsm.metadata['characteristics_ch1'][3].split(": ")[1])
    char["sample_name"].append(gsm.name)

print(char["stage"][3], char["time"][3], char["sex"][3], char["sample_name"][3])

GPL = gse.gpls.values()[0]
pivoted_samples = gse.pivot_samples('VALUE')
pivoted_samples.set_index(GPL.table.SPOT_ID, inplace=True)

#pivoted_samples.hist()

strata = pd.read_csv("phylostrata.txt", sep="\t", header=None)
strata.columns = ["GeneID", "ProbeID", "age"]
strata.set_index("ProbeID", inplace=True)

matched_data = pivoted_samples.join(strata, how="inner").groupby(level=0).last()

unique_data = matched_data.groupby("GeneID").mean()

char_pd = pd.DataFrame(char, index=char["sample_name"])
mixed = char_pd[char_pd.sex == "mixed"].sample_name.tolist()
mixed += char_pd[char_pd.sex == "female"].sample_name.tolist()

char_pd["timing_number"] = 0
time_stamps = char_pd.time.unique()
for i in xrange(len(time_stamps)):
    char_pd.loc[char_pd.time == time_stamps[i], "timing_number"] = i + 1

experiment_index = char_pd[char_pd.index.isin(mixed)].reset_index().groupby("timing_number")["index"].apply(lambda x: np.array(x))

set_mean = {}
stages = []
for d, col_list in experiment_index.iteritems():
    set_mean[d] = unique_data[col_list].mean(axis=1)
    stages.append(char_pd[char_pd.index.isin(col_list)].stage[0])

mean_data = pd.DataFrame(set_mean)
mean_data.plot()

# ===============
# Calculating TAI:

TAI = []
for s in mean_data.iteritems():
    expr_sum = sum(s[1])
    gcount = 0
    gene_TAI = 0
    for r in s[1]:
        gene_TAI += r*unique_data.age[gcount]
        gcount += 1
    TAI.append(gene_TAI/expr_sum)

plt.plot(TAI)
plt.xlabel()

# =============
# With log transformation
TAI = []
for s in np.log(mean_data).iteritems():
    expr_sum = sum(s[1])
    gcount = 0
    gene_TAI = 0
    for r in s[1]:
        gene_TAI += r*unique_data.age[gcount]
        gcount += 1
    TAI.append(gene_TAI/expr_sum)
plt.plot(TAI)

# =============

# ==============
# Plotting all data in a single histogram
all_expr = mean_data.values.tolist()
allflat = [item for sublist in all_expr for item in sublist]
plt.hist(allflat, bins=200)

# With log:

log_exp = np.log(mean_data.values)
log_expr = log_exp.tolist()
logflat = [item for sublist in log_exp for item in sublist]
plt.hist(logflat, bins=200)
# ==============
# Same thing with matrix multiplication:

age_indices = unique_data.age
expression_data = mean_data.values
product = np.dot(expression_data.T, age_indices)
mean_expression = expression_data.T.sum(1)
TAI = np.divide(product, mean_expression)

plt.plot(TAI)
plt.show()
plt.savefig("TAI.png")

TAI.tolist()