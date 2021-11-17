## Scripts dealing with VME analysis

### Data preparation

1. `merge_biigle_reports.py`: Merge BIIGLE reports by creating a single BIIGLE report for all survey and all labelling trees of interest.
2. `curate_biigle_report.py`: Cleans the BIIGLE report, removes rows with missing info (e.g., width, area, longitude). Saves a csv file where each row is an image, columns are: XX ... (sum of annotations present in the current image for this label).
3. `curate_coralnet_report.py`: Uses [ReadIn_Circumpolar_Annotation_Data_Prep](https://github.com/Southern-Ocean-Biodiversity-Mapping/ARC_Data/blob/main/ReadIn_Circumpolar_Annotation_Data_Prep.Rmd) csv file, selects VME morpho taxa.
4. `combine_coralnet_biigle.py`: Uses the output of `curate_biigle_report.py` and `curate_coralnet_report.py` and combines them. Accounts for different spatial coverage between BIIGLE and CoralNet (annotated images).
5. `csv_2_raster.R`: Converts the csv data, where abundance is aggregated, to raster, where abundance is aggregated per raster cell. Raster resolution can be changed. Saves a tif file where each layer contains the abundance data of each VME morpho-taxa, as well as metadata such as survey ID, imaged area.

#### TODO

4. `combine_coralnet_biigle.py`:
- [ ] Add infos: acquisition method, date --> waiting for Jan to propage info
- [ ] Add estimated area info --> waiting for Jan to propage info
- [ ] Check area_pix
5. `csv_2_raster.R`:
- [ ] Account for different spatial coverage between BIIGLE839 vs BIIGLE254
- [ ] Add missing covariates

### VME index

1. `run.R`: VME index mapping based on CCAMLR records. `config.R` contains the parameters needed for this analysis.
2. `ccamlr_vme_taxa_prevalence.py`: compare the `VMESpecimenWeight` and `VMESpecimenCount` for each VME taxon across the VME risk areas
![VMESpecimenCount_tot](https://user-images.githubusercontent.com/14353425/133349457-985fd82c-a1a8-4d67-9e72-a3299cbd73ac.png)

3. `generate_figure_1.R`: Show data distribution around Antarctica and sampling effort.
![image](https://user-images.githubusercontent.com/14353425/134629430-ecec814d-0b02-494e-a5ca-982b4c7c2895.png)
![image](https://user-images.githubusercontent.com/14353425/134629403-465f24e0-1e51-4201-a46c-7fa4cfd29457.png)

4. `generate_figure_2.R`: Looks at the variability in taxa abundance

![image](https://user-images.githubusercontent.com/14353425/135571117-8d0fabe0-23c9-4173-9dea-986707494bf7.png)
