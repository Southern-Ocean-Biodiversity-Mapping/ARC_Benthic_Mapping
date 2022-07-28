import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch

sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })

path_abd = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/20220523085637_020/df_vme_idx.csv"
path_div = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/20220523085350_020/df_vme_idx.csv"
path_conf = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_confidence_index.csv"
path_sea = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_sea_area.csv"
folder_out = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220614174710_021"

df_abd = pd.read_csv(path_abd)[["cellID", "VME index", "VME index category"]]
df_div = pd.read_csv(path_div)[["cellID", "VME index", "VME index category"]]
df_conf = pd.read_csv(path_conf)[["cellID", "confidence_index"]]
df_sea = pd.read_csv(path_sea)[["cellID", "section_name"]]
print(df_abd.head(5))
print(df_div.head(5))
print(df_conf.head(5))
print(df_sea.head(5))

print(df_div["VME index"].min(), df_div["VME index"].max())
print(df_div["VME index"].quantile(0.95), df_div["VME index"].quantile(0.8), df_div["VME index"].quantile(0.99))
div_perc95 = df_div["VME index"].quantile(0.95)
df_div_sea = pd.merge(df_div, df_sea, on="cellID")
print(df_div_sea[df_div_sea["VME index"] >= div_perc95].section_name.unique())

print(df_abd["VME index"].min(), df_abd["VME index"].max())
print(df_abd["VME index"].quantile(0.95), df_abd["VME index"].quantile(0.8), df_abd["VME index"].quantile(0.99))
abd_perc95 = df_abd["VME index"].quantile(0.95)
df_abd_sea = pd.merge(df_abd, df_sea, on="cellID")
print(df_abd_sea[df_abd_sea["VME index"] >= abd_perc95].section_name.unique())

#exit()

#clrs_vmc_abd_idx = []
#for v in df_abd["VME index category"].to_list():
#    if v == 1:
#        clrs_vmc_abd_idx.append("yellow")
#    elif v == 2:
#        clrs_vmc_abd_idx.append("orange")
#    else:
#        clrs_vmc_abd_idx.append("red")
#g = sns.barplot(data=df_abd, x="VME index", palette=clrs_vmc_abd_idx) #color="red")
sns.set_style(rc={'patch.force_edgecolor':True,
                   'patch.edgecolor': 'black'})

g = sns.displot(data=df_abd, x="VME index", hue="VME index category", height=5, aspect=2, palette=["yellow", "orange", "red"]) #color="red")
#g.set_xlabels("VMC abundance index")
g.savefig(os.path.join(folder_out, 'dist_abd.png'), dpi=300)
g = sns.displot(data=df_div, x="VME index", hue="VME index category", height=5, aspect=2, palette=["yellow", "orange", "red"])
#g.set_xlabels("VMC diversity index")
g.savefig(os.path.join(folder_out, 'dist_div.png'), dpi=300)
df_conf = df_conf[df_conf.confidence_index > 0]
def category_conf(x):
    if x <= 5:
        return 1
    elif x <= 7:
        return 2
    else:
        return 3
df_conf["conf category"] = df_conf.confidence_index.apply(lambda x: category_conf(x))
print(df_conf.sample(5))

g = sns.displot(data=df_conf, x="confidence_index", hue="conf category", height=5, aspect=2, palette=["black", "grey", "white"])
#g.set_xlabels("Confidence index")
g.savefig(os.path.join(folder_out, 'dist_conf.png'), dpi=300)