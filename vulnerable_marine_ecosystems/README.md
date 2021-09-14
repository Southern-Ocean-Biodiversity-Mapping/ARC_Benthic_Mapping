## Scripts dealing with VME analysis

Work in progress :-)

### Data preparation

1. `curate_biigle_report.py`: Cleans the BIIGLE report, removes rows with missing info (e.g., width, area, longitude). Saves a csv file where each row is an image, columns are: survey, image area, coordinates, and the coverage in m2 of each VME morpho-taxa (sum of annotations present in the current image for this label).
2. `curate_coralnet_report.py`: TODO
3. `combine_coralnet_biigle.py`: TODO
4. `biigle_2_raster.R`: Converts the csv data, where abundance is aggregated per image, to raster, where abundance is aggregated per raster cell. Raster resolution can be changed. Saves a tif file where each layer contains the abundance data of each VME morpho-taxa, as well as metadata such as survey ID, imaged area.

#### TODO

1. `curate_biigle_reports.py`
- [ ] Add infos: Image quality, acquisition method
- [ ] Fetch missing infos: Long, Lat, Width, Height, Area
- [ ] Account for empty images
2. `curate_coralnet_report.py`: TODO
3. `combine_coralnet_biigle.py`: TODO
4. `biigle_2_raster.R`:
- [ ] Rename to `csv_2_raster.R`

### VME index

TODO: show prelim results
