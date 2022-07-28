import os
import jenkspy
import pyreadr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch

sns.set_style("whitegrid", {
    'grid.linestyle': '--'
 })

path_bio_data = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_sea_area.csv"
path_metadata = "C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Data.Rdata"
path_out = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_confidence_index.csv"

folder_plot_out = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/20220330135843_009"

df_bio = pd.read_csv(path_bio_data)
df_im_metadata = pyreadr.read_r(path_metadata)["image_metadata"].reset_index()

df_im_metadata = df_im_metadata[df_im_metadata["counts"] == "yes"]

df_bio.drop(columns=["section_col", "section_name"], inplace=True)
df_im_metadata.drop(columns=["gear", "counts", "cover", "area_source", "survey", "proj_coord_x", "proj_coord_y",
                             "transectID", "lon", "lat", "Filename.standardised", "rownames"], inplace=True)

df_bio = df_bio.astype({"cellID": int})
df_im_metadata = df_im_metadata.astype({"cellID": int})

for idx, row in df_bio.iterrows():
    cell_id_cur = row["cellID"]
    df_im_metadata_cur = df_im_metadata[df_im_metadata["cellID"] == cell_id_cur]

    if df_im_metadata_cur.year.isnull().all():
        df_im_metadata_cur.loc[df_im_metadata_cur.index, "year"] = 2000

    if df_im_metadata_cur.area.isnull().all() or df_im_metadata_cur.image_quality_score.isnull().all():
        df_bio.loc[idx, "confidence_index"] = -1
    else:
        year_cur = df_im_metadata_cur["year"].mean()
        im_q_cur = df_im_metadata_cur["image_quality_score"].mean()
        area_cur = df_im_metadata_cur["area"].sum()
        if str(year_cur) == "nan":
            print(df_im_metadata_cur)
        df_bio.loc[idx, "age"] = 2022 - year_cur if str(year_cur) != "nan" else year_cur
        df_bio.loc[idx, "im_quality"] = im_q_cur
        df_bio.loc[idx, "sampled_portion"] = area_cur * 100. / (500 * 500)

print("Number of cells where it is not possible to compute the confidence index because missing data: {} ...".format(len(df_bio[df_bio["confidence_index"] == -1])))
print(df_bio[df_bio["confidence_index"] == -1])

n_breaks = 3

print("Breaking image quality...")
breaks_iq = jenkspy.jenks_breaks(df_bio["im_quality"],
                              nb_class=n_breaks)
df_bio["im_quality_score"] = pd.cut(df_bio["im_quality"],
                          bins=breaks_iq,
                          labels=[ll for ll in range(1, n_breaks+1)],
                          include_lowest=True)
df_bio["im_quality_score"].replace(to_replace=[ll for ll in range(1, n_breaks + 1)],
                         value=[ll for ll in range(1, n_breaks + 1)], inplace=True)
print(breaks_iq)
print(df_bio["im_quality_score"].value_counts())
print("")

print("Breaking sampled_portion...")
breaks_sa = jenkspy.jenks_breaks(df_bio["sampled_portion"],
                              nb_class=n_breaks)
df_bio["sampled_portion_score"] = pd.cut(df_bio["sampled_portion"],
                          bins=breaks_sa,
                          labels=[ll for ll in range(1, n_breaks+1)],
                          include_lowest=True)
df_bio["sampled_portion_score"].replace(to_replace=[ll for ll in range(1, n_breaks + 1)],
                         value=[ll for ll in range(1, n_breaks + 1)], inplace=True)
print(breaks_sa)
print(df_bio["sampled_portion_score"].value_counts())
print("")

print("Breaking age...")
breaks_a = jenkspy.jenks_breaks(df_bio["age"],
                              nb_class=n_breaks)
df_bio["age_score"] = pd.cut(df_bio["age"],
                          bins=breaks_a,
                          labels=[ll for ll in list(reversed(range(1, n_breaks + 1)))],
                          include_lowest=True)
df_bio["age_score"].replace(to_replace=[ll for ll in list(reversed(range(1, n_breaks + 1)))],
                         value=[ll for ll in list(reversed(range(1, n_breaks + 1)))], inplace=True)
print(breaks_a)
print(df_bio["age_score"].value_counts())
print("")

df_bio.loc[df_bio[df_bio["confidence_index"] != -1].index, "confidence_index"] = df_bio[["im_quality_score", "sampled_portion_score", "age_score"]].sum(axis=1)

print(df_bio.sample(10))

g = sns.displot(data=df_bio, x="sampled_portion", height=5, aspect=2, color="violet")
g.set_xlabels("Sampled area per cell [%]")
plt.axvline(x=breaks_sa[1], color='black', lw=2)
plt.axvline(x=breaks_sa[2], color='black', lw=2)
g.savefig(os.path.join(folder_plot_out, 'dist_conf_sampling.png'), dpi=300)

g = sns.displot(data=df_bio, x="age", height=5, aspect=2, color="mediumvioletred")
g.set_xlabels("Averaged imagery age per cell [%]")
plt.axvline(x=breaks_a[1], color='black', lw=2)
plt.axvline(x=breaks_a[2], color='black', lw=2)
g.savefig(os.path.join(folder_plot_out, 'dist_conf_age.png'), dpi=300)

g = sns.displot(data=df_bio, x="im_quality", height=5, aspect=2, color="magenta")
g.set_xlabels("Averaged imagery quality per cell [%]")
plt.axvline(x=breaks_iq[1], color='black', lw=2)
plt.axvline(x=breaks_iq[2], color='black', lw=2)
g.savefig(os.path.join(folder_plot_out, 'dist_conf_imQuality.png'), dpi=300)

exit()
df_bio.to_csv(path_out, index=False)
