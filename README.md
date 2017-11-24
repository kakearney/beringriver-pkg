
# Freshwater input to the Bering Sea, 1950-2017


Author: Kelly Kearney


This repository is provided as supplementary material to the NOAA technical report of the same title.


The primary functions in this repository are:



  - watersheds.m: Analyzes Watershed Boundary Dataset hydrologic units to   determine which ones are part of the Bering Sea/Gulf of Alaska   watersheds, and finds any USGS streamflow measurement stations   located within them
  - riversfromwatershed.m: Chooses which streamflow station measurements   to use for river mouth unit, and recontructs a single timeseries for   each.
  - riversToRunoff.m: Reformats freshwater discharge timeseries into a   spatially-resolved grid suitable for use as a ROMS forcing file.


## Contents

            
- Getting started        
- Useage

## Getting started


**Prerequisites**


This collection of functions requires Matlab R2015b or later, along with the Statistics and Machine Learning Toolbox, Mapping Toolbox, and Image Processing Toolbox.  The following additional code utilities are also used at various points in the code:



  - [ConsoleProgressBar](https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar): by Evgeny Pr, utility for text progress bar
  - [Shape Languge Modeling Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling): by John D'Errico, Least squares spline modeling using shape primitives
  - [borders](https://www.mathworks.com/matlabcentral/fileexchange/50390-borders): by Chad Greene, plots geographic borders
  - [consolidator](https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator): by John D'Errico, consolidates arrays
  - [distance2curve](https://www.mathworks.com/matlabcentral/fileexchange/34869-distance2curve): by John D'Errico, Find the closest point on a (n-dimensional) curve to any given point or set of points
  - [inpaintn](https://www.mathworks.com/matlabcentral/fileexchange/27994-inpaint-over-missing-data-in-1-d--2-d--3-d--nd-arrays): by Damien Garcia, inpaint over missing data * [JSONlab](https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files): by Qianqian Fang, encode/decode JSON files
  - [mergestruct](https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct): (renamed from catstruct.m) by Jos, merge fields from structure arrays * [plotSpread](https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points--beeswarm-plot-): by Jonas, plots a beeswarm plot
  - [subaxis](https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot): by Aslak Grinsted, better multi-axis layout
  - [textLoc](https://www.mathworks.com/matlabcentral/fileexchange/17151-textloc): by Ben Barrowes, intuitive text positioning
  - [roms_matlab](https://www.myroms.org): ROMS Matlab toolbox, available through the SVN repositories at myroms.org
  - BeringBESTNPZ: by Kelly Kearney, informal collection of tools related to the Bering 10K ROMS domain (contact me for this)
  - [aggregate](https://github.com/kakearney/aggregate-pkg): by Kelly Kearney, aggregates values in a matrix
  - [cellstr2](https://github.com/kakearney/cellstr2-pkg): by Kelly Kearney, extention of `cellstr` function
  - [cptcmap](https://github.com/kakearney/cptcmap-pkg): by Kelly Kearney, apply colormaps from color palette files
  - [dirfull](https://github.com/kakearney/dirfull-pkg): by Kelly Kearney, extension of `dir` function
  - [inpolygons](https://github.com/kakearney/inpolygons-pkg): by Kelly Kearney, extension of `inpolygon` function
  - [joinsegments](https://github.com/kakearney/joinsegments-pkg): by Kelly Kearney, connects NaN-delimited polygon/polyline segments
  - [mask2poly](https://github.com/kakearney/mask2poly-pkg): by Kelly Kearney, creates outline of logical mask
  - [minmax](https://github.com/kakearney/minmax-pkg): by Kelly Kearney, returns min and max of an array
  - [multitextloc](https://github.com/kakearney/multitextloc-pkg): by Kelly Kearney, extension of `textLoc` function
  - [ncreads](https://github.com/kakearney/ncreads-pkg): by Kelly Kearney, extension of `ncread` function
  - [offsetaxis](https://github.com/kakearney/offsetaxis-pkg): by Kelly Kearney, plots axis offset from plot area
  - [plotgrid](https://github.com/kakearney/plotgrid-pkg): by Kelly Kearney, sets up and/or plots to a grid of axes
  - [regexpfound](https://github.com/kakearney/regexpfound-pkg): by Kelly Kearney, return logical mask from `regexp` function
  - [shapeprjread](https://github.com/kakearney/shapeprjread-pkg): by Kelly Kearney, extension of `shaperead` function for projected shapefiles
  - [vlookup](https://github.com/kakearney/vlookup-pkg): by Kelly Kearney, implements lookup tables

In addition to these toolboxes, several large datasets are required that are not bundled here.  These include:



  - The Major Rivers data available through the Alaska State Geo-Spatial    Clearinghouse ([http://www.asgdc.state.ak.us/](http://www.asgdc.state.ak.us/)).  It includes selected    rivers from the USGS 1:2,000,000 Digital Line Graphs (DLG).
  - National Hydrography Dataset Hydrography Feature Dataset    (NHDFlowlines): 1:24000 scale flowlines, ([https://nhd.usgs.gov/](https://nhd.usgs.gov/))
  - Natural Earth 1:10m Rivers &amp; Lakes centerlines, ([http://www.naturalearthdata.com/](http://www.naturalearthdata.com/))

**Downloading and installation**


This code can be downloaded from [Github](https://github.com/kakearney/beringriver-pkg)


**Matlab Search Path**


This package relies on a number of external toolboxes (some my own, some third-party).  These must be added to your Matlab Search path (via `addpath`, `pathtool`, etc.):



## Useage


See comments in the three primary functions for a description of use. Note that this code is provided primarily for reference. For further information, or to aquire the exact datasets required   to run this code verbatim, please contact the author.



<sub>[Published with MATLAB R2016b]("http://www.mathworks.com/products/matlab/")</sub>