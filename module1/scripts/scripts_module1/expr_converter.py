import GEOparse
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

