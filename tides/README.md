# Tide Gauge Module

## Introduction

The Tide Gauge module was introduced in 

  Hawkins R., Bodin T., Sambridge M., Choblet G. and Husson L.,
  "Trans-dimensional surface reconstruction with different classes of parameterization",
  Geochemistry, Geophysics, Geosystems,
  2019

with a synthetic example around Tasmania. The synthetic examples appearing in the above
paper can be recreated using the Makefile within tas_synthetic sub-directory. The actual
results in the paper used a cluster for processing and the slurm scripts used can be
found in the tas_synthetic/transcale directory. These are cluster specific but may serve
as a guide.

## Virtual Tide Gauges

The Tide Gauge module has applied to real data (Tide Gauge, GPS and Satellite Radar
Altimetry) in

  Hawkins R., Husson L., Choblet G., Bodin T. and Pfeffer J.,
  "Virtual tide gauges for predicting relative sea level rise",
  JGR: Solid Earth,
  2019 (submitted)

### Data

The processed data as used in the above study appear in the jgr2019/data directory.

For each region there are four files:
  <region>_combined.txt
  <region>_regression_sea.txt
  <region>_regression_gps.txt
  <region>_regression_tides.txt

The format for the <region>_regression_*.txt files are simple text files where
the first line is the total number of points, then for each point there
are space seperated values for
<longitude> <latitude> <rate> <rate error>

The format for the <region>_combined.txt files are similar except for each point
there is a added "type" field, ie
<longitude> <latitude> <type> <rate> <rate error>

Which defines if the observation is a Tide Gauge (1), GPS (2) or Sea surface (3).

### Running

In the data_jgr2019/example directory there is a single processor inversion of
data for the Australian region. This is not representative of the final results
but is used to demonstrate the operation of the method on a small dataset.

The actual results from the paper are generated on a cluster and the
slurm submission scripts used are located in the data_jgr2019/slurm
subdirectory.

### Acknowledgements

These data files are constructed from the following data sets :

Tide Gauge Data from PSMSL
https://www.psmsl.org/

GPS Data from Nevada Geodetic Laboratory
http://geodesy.unr.edu/

Satellite Radar Altimetry from Copernicus Marine environment monitoring service
http://marine.copernicus.eu/

Pressure data from European Centre for Medium-Range Weather Forecasts
https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim



