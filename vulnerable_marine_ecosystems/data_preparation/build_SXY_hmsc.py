import os
import pyreadr
import pandas as pd

path_S = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_confidence_index.csv"
path_X = "C:/Users/cgros/code/IMAS/ARC_Data/annotation/Circumpolar_Annotation_Env_Data.RData"
path_Y = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/biodata_step4.csv"
path_out = "C:/Users/cgros/code/IMAS/ARC_Benthic_Mapping/vulnerable_marine_ecosystems/species_distribution_modelling/hmsc/data/SXY.csv"

df_S = pd.read_csv(path_S)
df_S_ = df_S[["cellID", "proj_coord_x", "proj_coord_y", "im_quality", "sampled_portion"]]
print(df_S_.head())
print(len(df_S_))

df_X = pyreadr.read_r(path_X)["cell_metadata_env"].reset_index()
df_X_ = df_X.drop(['index', 'lon', 'lat', 'proj_coord_x', 'proj_coord_y', 'cover_N', 'counts_N',
                   'cover_area', 'counts_area', 'cover_cells_survey', 'cover_cells_transect1', 'cover_cells_transect2',
                   'cover_cells_transect3', 'counts_cells_survey', 'counts_cells_transect1', 'counts_cells_transect2',
                   'counts_cells_transect3', 'image_quality_score', 'year', 'gear'], axis=1)
df_X_["cellID"] = pd.to_numeric(df_X_["cellID"])
print(df_X_.head())
print(df_X_.keys())
print(len(df_X_))

df_Y = pd.read_csv(path_Y)
df_Y_ = df_Y
print(df_Y_.head())
print(len(df_Y_))

df_SX = pd.merge(left=df_S_, right=df_X_, left_on="cellID", right_on="cellID")
df_SXY = pd.merge(left=df_SX, right=df_Y_, left_on="cellID", right_on="cellID")
print(df_SXY.head())
print(len(df_SXY))

df_SXY.to_csv(path_out, index=False)