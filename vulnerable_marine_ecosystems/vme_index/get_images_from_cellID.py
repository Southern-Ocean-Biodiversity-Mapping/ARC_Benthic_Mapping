import os
import pyreadr
import numpy as np
import pandas as pd

path_abd = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/20220523085637_020/df_vme_idx.csv"
path_div = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/vme_index/20220523085350_020/df_vme_idx.csv"
path_conf = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_confidence_index.csv"
path_sea = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_sea_area.csv"
path_im = "C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Data.Rdata"

df_abd = pd.read_csv(path_abd)[["cellID", "VME index"]]
df_div = pd.read_csv(path_div)[["cellID", "VME index"]]
df_conf = pd.read_csv(path_conf)
df_sea = pd.read_csv(path_sea)
df_im = pyreadr.read_r(path_im)["image_metadata"].reset_index()[["Filename.standardised", "cellID"]]
df_im["cellID"] = df_im["cellID"].astype(int)
df_conf["cellID"] = df_conf["cellID"].astype(int)

div_perc95 = df_div["VME index"].quantile(0.95)
div_perc50 = df_div["VME index"].quantile(0.5)
df_div_sea = pd.merge(df_div, df_sea, on="cellID")

#abd_perc95 = df_abd["VME index"].quantile(0.95)
#abd_perc50 = df_abd["VME index"].quantile(0.5)
#df_abd_sea = pd.merge(df_abd, df_sea, on="cellID")

#conf_perc95 = df_conf["confidence_index"].quantile(0.95)
#df_conf_sea = pd.merge(df_conf, df_sea, on="cellID")

#print("abd 95: {} abd 50: {}".format(abd_perc95, abd_perc50))
#print("div 95: {} div 50: {}".format(div_perc95, div_perc50))


def print_info(cellID_lst):
    lst_abd, lst_div, lst_conf = [], [], []
    for c in cellID_lst:
        print("\n" + str(c))
        abd_ = df_abd[df_abd["cellID"]==c]["VME index"].to_list()[0]
        div_ = df_div[df_div["cellID"]==c]["VME index"].to_list()[0]
        conf_ = df_conf[df_conf["cellID"] == c]["confidence_index"].to_list()[0]
        print("Div: {} \t Abd: {} \t Conf: {} \t Sea: {}".format(div_,
                                                                 abd_,
                                                                 conf_,
                                                                 df_sea[df_sea["cellID"]==c]["section_name"].to_list()[0]
                                                                 ))
        print("lon: {} \t lat: {}".format(round(df_sea[df_sea["cellID"]==c]["lon"].to_list()[0]),
                                          round(df_sea[df_sea["cellID"]==c]["lat"].to_list()[0])))
        print("age: {} \t im_quality {} \t sampled_portion {}".format(df_conf[df_conf["cellID"]==c]["age"].to_list()[0],
                                                                      df_conf[df_conf["cellID"]==c]["im_quality"].to_list()[0],
                                                                      df_conf[df_conf["cellID"]==c]["sampled_portion"].to_list()[0]))
        print("\n\t" + "\n\t".join(df_im[df_im.cellID == c]["Filename.standardised"].to_list()))

        lst_abd.append(abd_)
        lst_div.append(div_)
        lst_conf.append(conf_)

    print("Averaged Div {:0.2f} +/- {:0.2f}".format(np.mean(lst_div), np.std(lst_div)))
    print("Averaged Abd {:0.2f} +/- {:0.2f}".format(np.mean(lst_abd), np.std(lst_abd)))
    print("Averaged Conf {:0.2f} +/- {:0.2f}".format(np.mean(lst_conf), np.std(lst_conf)))

#FSR 881
#print_info([131969604])
#print_info([130849375])

#FSR 881 bis
#print_info([142625061])
#print_info([142451692])

#Peninsula
#print_info([46354318])
#print_info([45580853])

# Riisen
#print_info([40583904])
#print_info([39863856])

# Figure 9
print_info([54808593, 54861920])
#print_info([46300975])

# Inside VME
#print_info([46794111, 46807446, 145294701, 144081411, 144108104, 138730426, 138743761, 138743760, 141038096, 141251443, 141264778, 141678124, 141664789, 142451692, 142451693, 142451694, 133356496, 133356495, 133343159, 131969604, 131982938, 131662873, 131649538, 131636203])

# Inside FSR881 bis
#print_info([142451692, 142451693, 142451694])
# Outside FSR881 bis
#print_info([142518378, 142531714, 142545049, 142611724, 142625059, 142625060, 142625061, 142545070, 142558405, 142571740])
exit()
cellID_div_perc95_abd_perc50 = [c for c in df_div[df_div["VME index"] > div_perc95]["cellID"].to_list()
                              if c in df_abd[df_abd["VME index"] < abd_perc50]["cellID"].to_list()]
cellID_abd_div_perc95 = [c for c in df_div[df_div["VME index"] > div_perc95]["cellID"].to_list()
                              if c in df_abd[df_abd["VME index"] > abd_perc95]["cellID"].to_list()]
cellID_div_perc50_abd_perc95 = [c for c in df_div[df_div["VME index"] < div_perc50]["cellID"].to_list()
                              if c in df_abd[df_abd["VME index"] > abd_perc95]["cellID"].to_list()]

print_info(cellID_lst=cellID_abd_div_perc95)
#print_info(cellID_lst=cellID_div_perc95_abd_perc50)
#print_info(cellID_lst=cellID_div_perc50_abd_perc95)

