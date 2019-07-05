# TransTessellate2D
Trans-dimensional tesselation for 2D problems

This software is released in conjunction with the manuscript:

  Hawkins R., Bodin T., Sambridge M., Choblet G. and Husson L.,
  "Trans-dimensional surface reconstruction with different classes of parameterization",
  Geochemistry, Geophysics, Geosystems,
  2019

and provides a frame work for Trans-dimensional inference for 2D cartesian problems.

It is provided with the code for the three examples presented in the manuscript above
but allows extension to custom forward models and likelihoods.

See the documentation/manual.tex Latex file for a tutorial/manual on
how to run and customize the code.

# Compilation

This software is written in C++ and requires

- GNU g++ Version 6.x or greater
- GNU fortran Version 6.x or greater
- GNU Make version 4.x or greater
- GNU Scientific Library (GSL) version 1 or 2
- OpenMPI version 1.10

The source code with example scripts and data is called
TransTessellate2D.tar.gz and once the environment is properly
configued, extraction and compilation can proceed as follows:

```
> tar -xzf TransTessellate2D.tar.gz
> cd TransTessellate2D
> make 
```

## Virtual Tide Gauges

The Tide Gauge module of this software package has applied to real data (Tide Gauge, GPS and Satellite Radar
Altimetry) in

  Hawkins R., Husson L., Choblet G., Bodin T. and Pfeffer J.,
  "Virtual tide gauges for predicting relative sea level rise",
  JGR: Solid Earth,
  2019 (submitted)

For details, and data, see the README.md in the tides subdirectory.



