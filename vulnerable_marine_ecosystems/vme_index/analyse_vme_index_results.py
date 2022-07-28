import os
import numpy as np
import pandas as pd

path_geo_data = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_sea_area.csv"
path_abd_data = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220315100421_009/df_vme_idx.csv"
path_div_data = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220315090028_009/df_vme_idx.csv"

df_geo = pd.read_csv(path_geo_data)
df_abd = pd.read_csv(path_abd_data)
df_div = pd.read_csv(path_div_data)

df_abd = pd.merge(df_abd, df_geo[["cellID", "section_name", "cluster"]], on="cellID", how="right")
df_div = pd.merge(df_div, df_geo[["cellID", "section_name", "cluster"]], on="cellID", how="right")

df_abd_save = pd.merge(df_abd, df_geo, on="cellID", how="right")
df_abd_save.to_csv("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220315100421_009/df_vme_idx_geo.csv",
                   index=False)
df_div_save = pd.merge(df_div, df_geo, on="cellID", how="right")
df_div_save.to_csv("C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220315090028_009/df_vme_idx_geo.csv",
                   index=False)

for index_name, df_ in zip(["Abundance", "Diversity"], [df_abd, df_div]):
    print("\n{} VME index ...".format(index_name))
    perc_5 = np.percentile(df_["VME index"], 5)
    perc_95 = np.percentile(df_["VME index"], 95)

    df_perc_5 = df_[df_["VME index"] <= perc_5][["cellID", "section_name", "cluster"]]
    df_perc_95 = df_[df_["VME index"] >= perc_95][["cellID", "section_name", "cluster"]]

    print("\tMin - Max: {:.2f} - {:.2f} ...".format(min(df_["VME index"]), max(df_["VME index"])))

    print("\tPercentile 5% value: {:.2f} ...".format(perc_5))
    print("\tcomprising {} cells ...".format(len(df_perc_5)))
    print("\tacross {} clusters ...".format(len(df_perc_5["cluster"].unique())))
    print("\tin the sections: {} ...".format(df_perc_5["section_name"].unique()))

    print("\tPercentile 95% value: {:.2f} ...".format(perc_95))
    print("\tcomprising {} cells ...".format(len(df_perc_95)))
    print("\tacross {} clusters ...".format(len(df_perc_95["cluster"].unique())))
    print("\tin the sections: {} ...".format(df_perc_95["section_name"].unique()))

    df_cluster = df_.groupby(['cluster'])[["VME index"]].mean().sort_values(by="VME index")
    top_5 = df_cluster.iloc[-5:, ].index.to_list()
    print(df_cluster.iloc[-5:, ])
    print(df_[df_["cluster"].isin(top_5)]["section_name"].unique())

