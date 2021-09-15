## Scripts dealing with VME analysis

### Data preparation

1. `curate_biigle_report.py`: Cleans the BIIGLE report, removes rows with missing info (e.g., width, area, longitude). Saves a csv file where each row is an image, columns are: survey, image area, coordinates, and the coverage in m2 of each VME morpho-taxa (sum of annotations present in the current image for this label).
2. `curate_coralnet_report.py`: Uses [ReadIn_Circumpolar_Annotation_Data_Prep](https://github.com/Southern-Ocean-Biodiversity-Mapping/ARC_Data/blob/main/ReadIn_Circumpolar_Annotation_Data_Prep.Rmd) csv file, selects VME morpho taxa.
3. `combine_coralnet_biigle.py`: Uses the output of `curate_biigle_report.py` and `curate_coralnet_report.py` and combines them. Accounts for different spatial coverage between BIIGLE and CoralNet (annotated images).
4. `csv_2_raster.R`: Converts the csv data, where abundance is aggregated per image, to raster, where abundance is aggregated per raster cell. Raster resolution can be changed. Saves a tif file where each layer contains the abundance data of each VME morpho-taxa, as well as metadata such as survey ID, imaged area.

#### TODO

1. `curate_biigle_reports.py`
- [ ] Add infos: Image quality, acquisition method
- [ ] Fetch missing infos: Long, Lat, Width, Height, Area
- [ ] Account for empty images
- [ ] Add Victor and Jan label tree
2. `curate_coralnet_report.py`:
- [ ] Discuss UBS_Sp and UBS_B
- [ ] Rename columns to match BIIGLE convention
4. `combine_coralnet_biigle.py`:
- [ ] Account for different spatial coverage between BIIGLE and CoralNet (annotated images).
5. `csv_2_raster.R`:
- [ ] XX

### VME index

1. `run.R`: VME index mapping based on CCAMLR records. `config.R` contains the parameters needed for this analysis.
2. `ccamlr_vme_taxa_prevalence.py`: compare the `VMESpecimenWeight` and `VMESpecimenCount` for each VME taxon across the VME risk areas

![VMESpecimenCount_tot](https://user-images.githubusercontent.com/14353425/133349457-985fd82c-a1a8-4d67-9e72-a3299cbd73ac.png)


