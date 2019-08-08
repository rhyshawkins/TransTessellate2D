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
```
  <region>_combined.txt
  <region>_regression_sea.txt
  <region>_regression_gps.txt
  <region>_regression_tides.txt
```

The format for the &lt;region&gt;\_regression\_*.txt files are simple text files where
the first line is the total number of points, then for each point there
are space seperated values for
```
<longitude> <latitude> <rate> <rate error>
```

The format for the &lt;region&gt;\_combined.txt files are similar except for each point
there is a added "type" field, ie
```
<longitude> <latitude> <type> <rate> <rate error>
```

Which defines if the observation is a Tide Gauge (1), GPS (2) or Sea surface (3).

### Running

In the jgr2019/example directory there is a single processor inversion of
data for the Australian region. This is not representative of the final results
but is used to demonstrate the operation of the method on a small dataset.

The individual steps for running are programmed into the Makefile in the
jgr2019/example directory. The Makefile has comments in it to describe the
various stages of inversions, but can be simply run to obtain results. For example,
running

```
> cd jgr2019/example
> make
```

Will run the demonstration single chain inversion of the Australian dataset,
and produce 3 pdf plots of the mean of the ensemble created, e.g. mean_tide.pdf
can be compared to Figure 6 (g) in the paper. This inversion will take approximately
1 - 2 minutes on a modern desktop computer. The plotting routines require Numpy,
Matplotlib, and Basemap to be available.

### Slurm Scripts

The actual results from the paper are generated on a cluster and the
slurm submission scripts used are located in the jgr2019/slurm
subdirectory.

### Raw results data

Under the jgr2019/results directory are the raw data files containing the
mean, std. devation, and ensemble histogram files that can be used for
reanalysis. Each data file has a prefix where the meaning of the prefix
is

- "joint" all data type inversion
- "notg"  no tide gauge (GPS and SRA only)
- "tg"    tide gauge only inversion

Some simple python 3 plotting scripts are included in the jgr2019/scripts
directory and two bash scripts are included for recreating the figures from the paper.
For example:

```
> cd jgr2019/results
> bash plotmean.sh europe
> bash plotcoasthistogram.sh southamerica
```

Will first show the mean RSL for the Europe region (Figure 2(e)), and the second will show
the virtual tide gauge coastline for South America (Figure 5(b)). The 

The format of the results files are simple space separate text files for images
that can be simply loaded in python with the command "numpy.loadtxt(filename)". In
each region, an extents.txt file is included which gives the longitude/latitude
range. For the format of the histogram files, see the scripts/plot_coast_histogram.py
python code.

The plotting routines require Numpy, Matplotlib, and Basemap to be available.

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



