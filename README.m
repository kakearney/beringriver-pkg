%% Freshwater input to the Bering Sea, 1950-2017
% Author: Kelly Kearney
%
% This repository is provided as supplementary material to the NOAA
% technical report of the same title.
%
% The primary functions in this repository are:
%
% * watersheds.m: Analyzes Watershed Boundary Dataset hydrologic units to 
%   determine which ones are part of the Bering Sea/Gulf of Alaska 
%   watersheds, and finds any USGS streamflow measurement stations 
%   located within them
% * riversfromwatershed.m: Chooses which streamflow station measurements
%   to use for river mouth unit, and recontructs a single timeseries for 
%   each.
% * riversToRunoff.m: Reformats freshwater discharge timeseries into a
%   spatially-resolved grid suitable for use as a ROMS forcing file.
%   
%% Getting started
%
% *Prerequisites*
%
% This collection of functions requires Matlab R2015b or later, along with
% the Statistics and Machine Learning Toolbox, Mapping Toolbox, and Image
% Processing Toolbox.  The following additional code utilities are also
% used at various points in the code:
%
% * <https://www.mathworks.com/matlabcentral/fileexchange/30297-consoleprogressbar ConsoleProgressBar>: by Evgeny Pr, utility for text progress bar
% * <https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling Shape Languge Modeling Toolbox>: by John D'Errico, Least squares
% spline modeling using shape primitives
% * <https://www.mathworks.com/matlabcentral/fileexchange/50390-borders
% borders>: by Chad Greene, plots geographic borders
% * <https://www.mathworks.com/matlabcentral/fileexchange/8354-consolidator consolidator>: by John D'Errico, consolidates arrays
% * <https://www.mathworks.com/matlabcentral/fileexchange/34869-distance2curve distance2curve>: by John D'Errico, Find the closest point on a (n-dimensional) curve to any given point or set of points
% * <https://www.mathworks.com/matlabcentral/fileexchange/27994-inpaint-over-missing-data-in-1-d--2-d--3-d--nd-arrays inpaintn>: by Damien Garcia, inpaint over missing data
% *
% <https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files
% JSONlab>: by Qianqian Fang, encode/decode JSON files
% * <https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct mergestruct>: (renamed from catstruct.m) by Jos, merge fields from
% structure arrays 
% *
% <https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points--beeswarm-plot-
% plotSpread>: by Jonas, plots a beeswarm plot
% * <https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot subaxis>: by Aslak Grinsted, better multi-axis layout
% * <https://www.mathworks.com/matlabcentral/fileexchange/17151-textloc textLoc>: by Ben Barrowes, intuitive text positioning
% * <https://www.myroms.org roms_matlab>: ROMS Matlab toolbox, available
% through the SVN repositories at myroms.org
% * BeringBESTNPZ: by Kelly Kearney, informal collection of tools
% related to the Bering 10K ROMS domain (contact me for this)
% * <https://github.com/kakearney/aggregate-pkg aggregate>: by Kelly
% Kearney, aggregates values in a matrix
% * <https://github.com/kakearney/cellstr2-pkg cellstr2>: by Kelly Kearney,
% extention of |cellstr| function
% * <https://github.com/kakearney/cptcmap-pkg cptcmap>: by Kelly Kearney,
% apply colormaps from color palette files
% * <https://github.com/kakearney/dirfull-pkg dirfull>: by Kelly Kearney,
% extension of |dir| function
% * <https://github.com/kakearney/inpolygons-pkg inpolygons>: by Kelly Kearney, extension of |inpolygon| function
% * <https://github.com/kakearney/joinsegments-pkg joinsegments>: by Kelly Kearney, connects NaN-delimited polygon/polyline segments 
% * <https://github.com/kakearney/mask2poly-pkg mask2poly>: by Kelly
% Kearney, creates outline of logical mask 
% * <https://github.com/kakearney/minmax-pkg minmax>: by Kelly Kearney,
% returns min and max of an array
% * <https://github.com/kakearney/multitextloc-pkg multitextloc>: by Kelly
% Kearney, extension of |textLoc| function
% * <https://github.com/kakearney/ncreads-pkg ncreads>: by Kelly Kearney,
% extension of |ncread| function
% * <https://github.com/kakearney/offsetaxis-pkg offsetaxis>: by Kelly
% Kearney, plots axis offset from plot area
% * <https://github.com/kakearney/plotgrid-pkg plotgrid>: by Kelly Kearney,
% sets up and/or plots to a grid of axes
% * <https://github.com/kakearney/regexpfound-pkg regexpfound>: by Kelly
% Kearney, return logical mask from |regexp| function
% * <https://github.com/kakearney/shapeprjread-pkg shapeprjread>: by Kelly
% Kearney, extension of |shaperead| function for projected shapefiles
% * <https://github.com/kakearney/vlookup-pkg vlookup>: by Kelly Kearney,
% implements lookup tables
%
% In addition to these toolboxes, several large datasets are required that
% are not bundled here.  These include:
%
% * The Major Rivers data available through the Alaska State Geo-Spatial
%    Clearinghouse (<http://www.asgdc.state.ak.us/>).  It includes selected
%    rivers from the USGS 1:2,000,000 Digital Line Graphs (DLG).
% * National Hydrography Dataset Hydrography Feature Dataset
%    (NHDFlowlines): 1:24000 scale flowlines, (<https://nhd.usgs.gov/>)
% * Natural Earth 1:10m Rivers & Lakes centerlines,
% (<http://www.naturalearthdata.com/>)
%
% *Downloading and installation*
%
% This code can be downloaded from <https://github.com/kakearney/beringriver-pkg Github>
%  
% *Matlab Search Path*
%
% This package relies on a number of external toolboxes (some my own, some
% third-party).  These must be added to your Matlab Search path (via
% |addpath|, |pathtool|, etc.):  
%
%% Useage
%
% See comments in the three primary functions for a description of use.
% Note that this code is provided primarily for reference. For further
% information, or to aquire the exact datasets required   to run this code
% verbatim, please contact the author.  




