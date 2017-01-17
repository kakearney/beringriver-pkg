%% Freshwater input to the Bering Sea, 1950--2016
% Author: Kelly Kearney
%
% This repository includes the code for the |buildberingrunoff.m| Matlab
% function, along with all dependent functions required to run it.  
%
% This function collects river runoff data (primarily from USGS) and uses
% it to reconstruct historical freshwater input to the Bering Sea, then
% builds an appropriate runoff forcing file for use with a ROMS model
% domain located in this region.
%
% This repository is provided as supplementary material to the NOAA
% technical report of the same title.
%
%% Getting started
%
% *Prerequisites*
%
% This function requires Matlab R2015b or later, along with the Statistics
% and Machine Learning Toolbox and Mapping Toolbox.  It also requires
% access to the <http://nco.sourceforge.net/ netCDF
% operators (NCO)> command line tools.
%
% *Downloading and installation*
%
% This code can be downloaded from <https://github.com/kakearney/beringriver-pkg Github>
%  
% *Matlab Search Path*
%
% This package relies on a number of external toolboxes (some my own, some
% third-party).  These have all been included in this repository. The
% following folders need to be added to your Matlab Search path (via
% |addpath|, |pathtool|, etc.):  
%
%  beringriver-pkg/FEX-function_handle
%  beringriver-pkg/SLMtools/SLMtools
%  beringriver-pkg/aggregate
%  beringriver-pkg/beringriver
%  beringriver-pkg/borders
%  beringriver-pkg/bufferm2
%  beringriver-pkg/cellstr2
%  beringriver-pkg/ginput2
%  beringriver-pkg/jsonlab-1.0
%  beringriver-pkg/mergestruct
%  beringriver-pkg/minmax
%  beringriver-pkg/ncreads
%  beringriver-pkg/regexpfound
%  beringriver-pkg/rgb_v3
%  beringriver-pkg/roms_matlab
%  beringriver-pkg/shapeprjread
%  beringriver-pkg/subgrid
%  beringriver-pkg/vlookup

%% Syntax

help buildberingrunoff


