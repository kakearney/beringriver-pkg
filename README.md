
# Freshwater input to the Bering Sea, 1950--2016


Author: Kelly Kearney


This repository includes the code for the `buildberingrunoff.m` Matlab function, along with all dependent functions required to run it.


This function collects river runoff data (primarily from USGS) and uses it to reconstruct historical freshwater input to the Bering Sea, then builds an appropriate runoff forcing file for use with a ROMS model domain located in this region.


This repository is provided as supplementary material to the NOAA technical report of the same title.



## Contents

            
- Getting started        
- Syntax

## Getting started


**Prerequisites**


This function requires Matlab R2015b or later, along with the Statistics and Machine Learning Toolbox and Mapping Toolbox.  It also requires access to the [netCDF operators (NCO)](http://nco.sourceforge.net/) command line tools.


**Downloading and installation**


This code can be downloaded from [Github](https://github.com/kakearney/beringriver-pkg)


**Matlab Search Path**


This package relies on a number of external toolboxes (some my own, some third-party).  These have all been included in this repository. The following folders need to be added to your Matlab Search path (via `addpath`, `pathtool`, etc.):



```
beringriver-pkg/FEX-function_handle
beringriver-pkg/SLMtools/SLMtools
beringriver-pkg/aggregate
beringriver-pkg/beringriver
beringriver-pkg/borders
beringriver-pkg/bufferm2
beringriver-pkg/cellstr2
beringriver-pkg/ginput2
beringriver-pkg/jsonlab-1.0
beringriver-pkg/mergestruct
beringriver-pkg/minmax
beringriver-pkg/ncreads
beringriver-pkg/regexpfound
beringriver-pkg/rgb_v3
beringriver-pkg/roms_matlab
beringriver-pkg/shapeprjread
beringriver-pkg/subgrid
beringriver-pkg/vlookup
```



## Syntax



```matlab
help buildberingrunoff
```




```
 BUILDBERINGRUNOFF  Run the river timeseries reconstruction process.
 
  [RivList, S, D, discharge, Dts, RusInv, precip] = ...
     buildberingrunoff(grdfile, outbase, varargin) 
 
  This function collects river runoff data (primarily from USGS) and uses
  it to reconstruct historical freshwater input to the Bering Sea, then
  builds an appropriate runoff forcing file for use with a ROMS model
  domain located in this region.
 
  While this function has been written to be easy to rerun periodically (as
  new data becomes available), there are several places in the code that
  were hard-coded based on a lot of manual analysis (especially in the
  "Build timeseries" code cell).  I recommend double-checking the
  assumptions in that section whenever rerunning... as stations are added,
  removed, and extended, the old fit routines may no longer be the best
  choice.
 
  This function loads several large external files (river coordinates from
  several sources, the Russian river dataset).  The paths to these files
  are set using optional input parameters whose defaults are currently set
  to match my computer; update these as necessary to run this function
  elsewhere.
 
  Input variables:
 
    grdfile:    full path to ROMS grid file
 
    outbase:    base name for output netcdf file.  Output file will be
                named [outbase].[yr1]-[yr2].efol[efol]km.updated[yymmdd],
                where yr1 and yr2 are the first and last year,
                respectively,  of river data included in the file, efol is
                the e-folding scale, and yymmdd is the date this function
                is run.
 
  Optional input variables (passed as parameter/value pairs):
 
    efol:       e-folding length scale (km) used to distribute runoff
                across grid cells [20]
 
    readcodes:  logical scalar indicating whether to download and read a
                new list of parameter codes from the USGS site.  These
                details are only occasionally changed by USGS, and I doubt
                they would ever change the code associated with discharge,
                but you may want to periodically refetch this data, just in
                case.  Codes are saved to the current folder using the name
                parameter_cd_query_XXX.txt, where XXX is one of the code
                categories.  If false, existing parameter_cd_query*.txt
                files in the current directory will be used. [false]
 
    readsites:  logical scalar indicating whether to query the USGS
                database for a list of all Alaskan monitoring sites (active
                and inactive) and the parameter codes measured at each.
                Data is saved to a file in the current directory named
                akriverlist.xml.  This step should be rerun any time you
                intend to download new discharge data, to check for any new
                sites becoming active.  If false, the exisiting
                akriverlist.xml file in the current directory will be used.
                [false]
 
    chooseriv:  logical scalar indicating whether to open the
                river-choosing window in order to set up the river names
                and mouth locations.  See chooserivers.m for more details.
                This step only needs to be rerun if you intend to add a new
                river to the dataset. [false]
 
    fetchdata:  logical scalar indicating whether to redownload discharge
                data from the USGS website.  If true, data will be queried
                and returned as JSON files, saved to a folder in the
                current directory named rivertimeseriesYYYYMMDD, where
                YYYYMMDD is today's date.  This step should be rerun
                whenever you want to add more recent data to the final
                forcing file. [false]
 
    choosedata: logical scalar indicating whether to redo the interactive
                data analysis.  This consists of plotting the timeseries
                associated with each river, and taking notes on how you
                want to use the data in the "Build timeseries" section of
                code.  This step should be repeated any time new data is
                fetched.  After this step, you may need to go in and change
                the hardcoded "Build timeseries" portion of this code.
                [false] 
 
    writefile:  logical scalar indicating whether to create a new netcdf
                forcing file holding the results of these calculations.
                Can be set to false if you just want to retrieve
                intermediate values (see output) without creating a new
                file. [false]
 
    nhdfile:    path to .mat file holding National Hydrography Database
                flowline data (see nhddata.m).  Only needed if chooseriv =
                true.
                ['nhddata.mat']
  
    majrivfile: path to shapefile of major Alaskan Rivers, downloaded from
                the Alaska State Geospatial
                Clearinghouse:http://www.asgdc.state.ak.us/ 
                ['/Volumes/Storage/AKGeospatialClearinghouse/export_1481830586108/mv_major_river_ln.shp']
 
    fsuinvfile: path to the list of rivers associated with the Russian
                river dataset (UCAR 553.2)
                ['/Volumes/Storage/UCAR_553p2_russiaRivers/fsu.inv']
 
    fsudatafile:path to the river flow data associated with the Russian
                river dataset 
                ['/Volumes/Storage/UCAR_553p2_russiaRivers/fsu2.txt']
 
    nerivfile:  path to the Natural Earth river and lake centerline
                shapefile (large scale 1:10m version)
                ['/Volumes/Storage/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp']
 
  Output variables:
 
    RivList:    table listing all rivers, their river mouth locations, and
                the USGS sites located along them
 
    S:          table listing all USGS sites (by name and site code), their
                location (lat, lon, and hydrographic unit code), and the
                dates over which each data measurement type is available at
                each station.
 
    D:          table of sites along the rivers selected for this study,
                including site name, site code, lat, lon, and the available
                time vs discharge timeseries for each site.
 
    tall:       ntime x 1 array, union of time values in all timeseries
                tables in D (sorted, but not evenly spaced)
 
    discharge:  ntime x nsite array, discharge data from D mapped onto an
                even time grid, with NaNs wherever data not available.
 
    Dts:        structure holding filled data:
 
                filled:         ntime x nriv array, where the columns
                                correspond to the rows of RivList 
                filltype:       ntime x nriv x 2, filltype(:,:,1) = method
                                used to fill in values in filled : 1 =
                                direct copy/paste, 2 = SLM-fitting, 0 = not
                                filled   filltype(:,:,2) = index of site
                                used   
                tmonth:         nmon x 1 datetime array, dates
                                corresponding to the rows of the 'month'
                                field  
                month:          nmon x nriv array, monthly averages of
                                'filled' field 
                clima:          12 x nriv array, climatological monthly
                                averages of 'filled' field 
                filltypemonth:  nmon x nriv array, primary method used in
                                each month-block 
                fillsitemonth:  nmon x nriv array, primary site used in
                                each month-block 
 
    RusInv:     table of russian river gauge details (data from fsuinv
                file, see above)  
 
    precip:     nxi x neta x ntime array of surface freshwater flux (i.e.
                runoff-as-precipitation), in kg/m^2/s.  This is the array
                saved as the Runoff variable in the output netCDF file. 


```



<sub>[Published with MATLAB R2016a]("http://www.mathworks.com/products/matlab/")</sub>