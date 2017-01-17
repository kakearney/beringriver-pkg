function [RivList, S, D, discharge, Dts] =  buildberingrunoff(grdfile, outbase, varargin)
%BUILDBERINGRUNOFF  Run the river timeseries reconstruction process.
%
% This function collects river runoff data (primarily from USGS) and uses
% it to reconstruct historical freshwater input to the Bering Sea, then
% builds an appropriate runoff forcing file for use with a ROMS model
% domain located in this region.
%
% While this function has been written to be easy to rerun periodically (as
% new data becomes available), there are several places in the code that
% were hard-coded based on a lot of manual analysis (especially in the
% "Build timeseries" code cell).  I recommend double-checking the
% assumptions in that section whenever rerunning... as stations are added,
% removed, and extended, the old fit routines may no longer be the best
% choice.
%
% This function loads several large external files (river coordinates from
% several sources, the Russian river dataset).  The paths to these files
% are set using optional input parameters whose defaults are currently set
% to match my computer; update these as necessary to run this function
% elsewhere.
%
% Input variables:
%
%   grdfile:    full path to ROMS grid file
%
%   outbase:    base name for output netcdf file.  Output file will be
%               named [outbase].[yr1]-[yr2].efol[efol]km.updated[yymmdd],
%               where yr1 and yr2 are the first and last year,
%               respectively,  of river data included in the file, efol is
%               the e-folding scale, and yymmdd is the date this function
%               is run.
%
% Optional input variables (passed as parameter/value pairs):
%
%   efol:       e-folding length scale (km) used to distribute runoff
%               across grid cells [20]
%
%   readcodes:  logical scalar indicating whether to download and read a
%               new list of parameter codes from the USGS site.  These
%               details are only occasionally changed by USGS, and I doubt
%               they would ever change the code associated with discharge,
%               but you may want to periodically refetch this data, just in
%               case.  Codes are saved to the current folder using the name
%               parameter_cd_query_XXX.txt, where XXX is one of the code
%               categories.  If false, existing parameter_cd_query*.txt
%               files in the current directory will be used. [false]
%
%   readsites:  logical scalar indicating whether to query the USGS
%               database for a list of all Alaskan monitoring sites (active
%               and inactive) and the parameter codes measured at each.
%               Data is saved to a file in the current directory named
%               akriverlist.xml.  This step should be rerun any time you
%               intend to download new discharge data, to check for any new
%               sites becoming active.  If false, the exisiting
%               akriverlist.xml file in the current directory will be used.
%               [false]
%
%   chooseriv:  logical scalar indicating whether to open the
%               river-choosing window in order to set up the river names
%               and mouth locations.  See chooserivers.m for more details.
%               This step only needs to be rerun if you intend to add a new
%               river to the dataset. [false]
%
%   fetchdata:  logical scalar indicating whether to redownload discharge
%               data from the USGS website.  If true, data will be queried
%               and returned as JSON files, saved to a folder in the
%               current directory named rivertimeseriesYYYYMMDD, where
%               YYYYMMDD is today's date.  This step should be rerun
%               whenever you want to add more recent data to the final
%               forcing file. [false]
%
%   choosedata: logical scalar indicating whether to redo the interactive
%               data analysis.  This consists of plotting the timeseries
%               associated with each river, and taking notes on how you
%               want to use the data in the "Build timeseries" section of
%               code.  This step should be repeated any time new data is
%               fetched.  After this step, you may need to go in and change
%               the hardcoded "Build timeseries" portion of this code.
%               [false] 
%
%   writefile:  logical scalar indicating whether to create a new netcdf
%               forcing file holding the results of these calculations.
%               Can be set to false if you just want to retrieve
%               intermediate values (see output) without creating a new
%               file. [false]
%
%   nhdfile:    path to .mat file holding National Hydrography Database
%               flowline data (see nhddata.m).  Only needed if chooseriv =
%               true.
%               ['nhddata.mat']
% 
%   majrivfile: path to shapefile of major Alaskan Rivers, downloaded from
%               the Alaska State Geospatial
%               Clearinghouse:http://www.asgdc.state.ak.us/ 
%               ['/Volumes/Storage/AKGeospatialClearinghouse/export_1481830586108/mv_major_river_ln.shp']
%
%   fsuinvfile: path to the list of rivers associated with the Russian
%               river dataset (UCAR 553.2)
%               ['/Volumes/Storage/UCAR_553p2_russiaRivers/fsu.inv']
%
%   fsudatafile:path to the river flow data associated with the Russian
%               river dataset 
%               ['/Volumes/Storage/UCAR_553p2_russiaRivers/fsu2.txt']
%
%   nerivfile:  path to the Natural Earth river and lake centerline
%               shapefile (large scale 1:10m version)
%               ['/Volumes/Storage/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp']
%
% Output variables:
%
%   RivList:    table listing all rivers, their river mouth locations, and
%               the USGS sites located along them
%
%   S:          table listing all USGS sites (by name and site code), their
%               location (lat, lon, and hydrographic unit code), and the
%               dates over which each data measurement type is available at
%               each station.
%
%   D:          table of sites along the rivers selected for this study,
%               including site name, site code, lat, lon, and the available
%               time vs discharge timeseries for each site.
%
%   tall:       ntime x 1 array, union of time values in all timeseries
%               tables in D (sorted, but not evenly spaced)
%
%   discharge:  ntime x nsite array, discharge data from D mapped onto an
%               even time grid, with NaNs wherever data not available.
%
%   Dts:        structure holding filled data:
%
%               filled:     ntime x nriv array, where the columns
%                           correspond to the rows of RivList 
%               filltype:   ntime x nriv x 2, 
%                           filltype(:,:,1) = method used to fill in values
%                           in filled : 1 =  direct copy/paste, 2 =
%                           SLM-fitting, 0 = not filled  
%                           filltype(:,:,2) = index of site used
%               month:      monthly averages of filled
%               clima:      climatological monthly averages of filled
%               filltypemonth:  primary method used in each month-block
%               fillsitemonth:  primary site used in each month-block

% Copyright 2016-2017 Kelly Kearney

%% User-set parameters
%
% outbase:  name to use as base for output netcdf file
% grdfile:  path to ROMS grid file
% efol:     e-folding scale used to distribute streamflow to grid cells

validateattributes(grdfile, {'char'}, {}, 'beringfreshwater', 'grdfile');
validateattributes(outbase, {'char'}, {}, 'beringfreshwater', 'outbase');

if ~exist(grdfile, 'file')
    error('Grid file (%s) not found');
end

p = inputParser;
p.addParameter('efol',        20,    @(x) validateattributes(x, {'numeric'},{'scalar'}));
p.addParameter('readcodes',   false, @(x) validateattributes(x, {'logical'},{'scalar'}));
p.addParameter('readsites',   false, @(x) validateattributes(x, {'logical'},{'scalar'}));
p.addParameter('chooseriv',   false, @(x) validateattributes(x, {'logical'},{'scalar'}));
p.addParameter('fetchdata',   false, @(x) validateattributes(x, {'logical'},{'scalar'}));
p.addParameter('choosedata',  false, @(x) validateattributes(x, {'logical'},{'scalar'}));
p.addParameter('writefile',   false, @(x) validateattributes(x, {'logical'},{'scalar'}));

p.addParameter('nhdfile',     'nhddata.mat', @(x) validateattributes(x, {'char'},{}));
p.addParameter('majrivfile',  '/Volumes/Storage/AKGeospatialClearinghouse/export_1481830586108/mv_major_river_ln.shp', @(x) validateattributes(x, {'char'},{}));
p.addParameter('fsuinvfile',  '/Volumes/Storage/UCAR_553p2_russiaRivers/fsu.inv', @(x) validateattributes(x, {'char'},{}));
p.addParameter('fsudatafile', '/Volumes/Storage/UCAR_553p2_russiaRivers/fsu2.txt', @(x) validateattributes(x, {'char'},{}));
p.addParameter('nerivfile',   '/Volumes/Storage/NaturalEarth/ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp', @(x) validateattributes(x, {'char'},{}));

p.parse(varargin{:});

Opt = p.Results;

% % Main ROMS Bering 10K grid
% 
% fol = '/Volumes/Storage/BeringROMS';
% grdfile = fullfile(fol, 'Bering_grid_withFeast.nc');
% outbase = 'riverrunoff-usgs-fsu';
% 
% % Al's Bristol Bay grid
% 
% % grdfile = 'roms_grd_drag.nc';
% % outbase = 'riverrunoff-usgs-fsu_grd_drag';
% 
% % E-folding scale (km)
% 
% efol = 20;
% 
% % Which actions should we redo?
% 
% readcodes = false; % true: query USGS database and save to file, false: use existing files
% readsites = false; % true: query USGS database and save to file, false: use existing files
% chooseriv = false; % true: use map to choose rivers, false: load from file
% fetchdata = false; % true: download and parse timeseries from NWIS database, false: load from file
% choosedata = false; % true: look at timeseries to choose best station, which ones need SLM fitting, etc, false: assume all is good to go
% writetofile = false; % true: write to file, false: stop after timeseries are constructed

%% Gather USGS river stations
%
% This cell creates 2 tables:
%
% Cd:   table of all USGS parameter codes, with descriptions, units, etc.
% S:    table of all stations, including some identifying info as well
%       as the most recent date (if any) for which any given parameter was
%       measured at that station

fprintf('Gathering USGS river stations...\n');

% Query for tables of parameter codes, save to file

cdtypes = {'PHY', 'BIO', 'INM', 'SED'};
cdfiles = cellfun(@(x) sprintf('parameter_cd_query_%s.txt', x), cdtypes, 'uni', 0);
ncdfile = length(cdfiles);

if Opt.readcodes
    for ii = 1:ncdfile
        cmd = sprintf('https://help.waterdata.usgs.gov/code/parameter_cd_query?fmt=rdb&group_cd=%s', cdtypes{ii});
        [f, status] = urlwrite(cmd, cdfiles{ii});
        if ~status
            error('Failed to read parameter codes: %s', cdtypes{ii});
        end
    end
end

cds = cell(ncdfile,1);
for ii = 1:ncdfile
    cds{ii} = readtable(cdfiles{ii}, 'delimiter', '\t', 'headerlines', 7);
    cds{ii} = cds{ii}(2:end,:);
end
Cd = cat(1, cds{:});

% Query for site list: daily values datasets in Alaska, streams only, with
% mean data for most recent date collected, all possible parameters

urlcmd = 'https://waterservices.usgs.gov/nwis/dv/?format=waterml,1.1&stateCd=ak&siteType=ST';
sitefile = 'akriverlist.xml';
if Opt.readsites
    [f, status] = urlwrite(urlcmd, sitefile);
    if ~status
        error('Failed to read list of sites');
    end
end

% Prepare to parse the site list xml file

dom = xmlread(sitefile);

% Which parameters are measured at each station?

tmp = dom.getElementsByTagName('ns1:variableCode');
nts = tmp.getLength;
vcode = cell(nts,1);
for ii = 1:nts
    vcode{ii} = char(tmp.item(ii-1).item(0).getData);
end

unqcd = unique(vcode);
ncd = length(unqcd);

opt = dom.getElementsByTagName('ns1:option');
ocode = cell(nts,1);
for ii = 1:nts
    ocode{ii} = char(opt.item(ii-1).item(0).getData);
end

% Read all data into a structure

S = struct;

idx = find(strcmp(ocode, 'Mean')) - 1;
n = length(idx);
[S.name, S.code, S.date, S.huc] = deal(cell(n,1));
[S.lat, S.lon] = deal(nan(n,1));

for ii = 1:n
    sinfo = tmp.item(idx(ii)).getParentNode.getParentNode;

    S.name{ii} = char(sinfo.getElementsByTagName('ns1:siteName').item(0).item(0).getData);
    S.code{ii} = char(sinfo.getElementsByTagName('ns1:siteCode').item(0).item(0).getData);
    S.lat(ii) = str2double(char(sinfo.getElementsByTagName('ns1:latitude').item(0).item(0).getData));
    S.lon(ii) = str2double(char(sinfo.getElementsByTagName('ns1:longitude').item(0).item(0).getData));
    
    huc = sinfo.getElementsByTagName('ns1:siteProperty');
    for ih = 1:huc.getLength
        propname = char(huc.item(ih-1).getAttribute('name'));
        if strcmp(propname, 'hucCd')
            S.huc{ii} = char(huc.item(ih-1).item(0).getData);
            break
        end
    end
    
    val = sinfo.getElementsByTagName('ns1:value');
    S.date{ii} = char(val.item(0).getAttribute('dateTime'));
end

% Group by station.  cd_XXXXX fields list most recent date that parameter
% is available at that station.

[stcd, is] = unique(S.code);
[~, stationidx] = ismember(S.code, stcd);
[~, vidx] = ismember(vcode(idx+1), unqcd);

tmp = cell(length(stcd), length(unqcd));
ind = sub2ind(size(tmp), stationidx, vidx);
tmp(ind) = S.date;

S = rmfield(S, 'date');
S = struct2table(S);
S = S(is,:);
cdcol = cellfun(@(x) sprintf('cd_%s', x), unqcd, 'uni', 0);
S = [S cell2table(tmp, 'VariableNames', cdcol)];

% Parse river/creek names

sname = regexp(S.name, '\s(AT|NR|BL|AB|ON|IN|[\d\.]*MI|[\d\.]* MI)\s', 'split');
sname = cellfun(@(x) x{1}, sname, 'uni', 0);

hasdis = ~cellfun('isempty', S.cd_00060);
islate = false(size(S,1),1);
islate(hasdis) = datenum(S.cd_00060(hasdis), 'yyyy-mm-ddTHH:MM:SS.FFF') > datenum(1950,1,1);

huc = unique(S.huc);
[~, hucidx] = ismember(S.huc, huc);

%% Identify rivers that flow into our ROMS domain
%
% This cell creates the RivList table, with names and river mouth
% coordinates of relevant rivers, as well as a list of all USGS stations on
% each river.
%
% I browsed through a *lot* of different river datasets looking for one
% that showed path coordinates for all the potential rivers monitored by
% the USGS gauges.  Most datasets (e.g. Natural Earth) only include the
% biggest ones.  In the end, I settled on two datasets:
% 1) The Major Rivers data available through the Alaska State Geo-Spatial
%    Clearinghouse (http://www.asgdc.state.ak.us/).  It includes selected
%    rivers from the USGS 1:2,000,000 Digital Line Graphs (DLG).
% 2) National Hydrography Dataset Hydrography Feature Dataset
%    (NHDFlowlines): 1:24000 scale flowlines, includes every little stream,
%    grouped by hydrographic unit. See nhddata.m for details of how I
%    preprocessed this data. 

fprintf('Gathering list of rivers...\n');

% Load river datasets

R = shapeprjread(Opt.majrivfile);

% Box around grid

Grd = ncreads(grdfile);

[~, grdfileshort, ex] = fileparts(grdfile);
grdfileshort = [grdfileshort ex];

boxfun = @(x) [x(1,:) x(:,end)' x(end,end:-1:1) x(end:-1:1,1)'];
boxlon = boxfun(wrapTo360(Grd.lon_rho));
boxlat = boxfun(Grd.lat_rho);
[boxlat, boxlon] = reducem(boxlat', boxlon');

hfig = figure; % buffer in projected coordinates to avoid dateline problems
usamap('alaska');
[x,y] = mfwdtran(boxlat, boxlon);
[x,y] = bufferm2('xy', x, y, 77500*2, 'out');
[boxlat2, boxlon2] = minvtran(x,y);
close(hfig);
boxlon2 = wrapTo360(boxlon2);

% Choose rivers

if Opt.chooseriv
    Nhd = load(Opt.nhdfile);
    chooserivers('beringrivermouths', S, R, Nhd, boxlat, boxlon, sname);
end

M = load('beringrivermouths');

% Final river list

RivList.name = upper(M.mname)';
RivList.lat = M.mlat';
RivList.lon = M.mlon';

nriv = length(RivList.lat);

% Match up sites with rivers

[tmp, sites] = aggregate(sname, S.code);

[tf, loc] = ismember(RivList.name, tmp);

rivsites = cell(nriv,1);
rivsites(tf) = sites(loc(tf));

% There are two Snake Rivers (of course there are...).  Sort those out.

issnake = strcmp(tmp, 'SNAKE R');
snakelat = vlookup(S, sites{issnake}, 'lat', 'code');
snakelon = vlookup(S, sites{issnake}, 'lon', 'code');

ridx = find(regexpfound(RivList.name, 'SNAKE'));

RivList.name(ridx) = cellstr(num2str((1:length(ridx))', 'SNAKE R (%d)'));
[rivsites{ridx}] = deal(cell(0));

for ii = 1:length(snakelat)
    d = distance(snakelat(ii), snakelon(ii), RivList.lat(ridx), RivList.lon(ridx));
    [~,imin] = min(d);
    rivsites{ridx(imin)} = [rivsites{ridx(imin)}; sites{issnake}(ii)];
end

% There also seems to be a second Wood River far inland

iswood = strcmp(tmp, 'WOOD R');
wtmp = sites{iswood};
woodlat = vlookup(S, wtmp, 'lat', 'code');
woodlon = vlookup(S, wtmp, 'lon', 'code');

ridx = find(regexpfound(RivList.name, 'WOOD'));
ell = referenceEllipsoid('earth');
d = distance(woodlat, woodlon, RivList.lat(ridx), RivList.lon(ridx), ell)/1000;
isclose = d < 100;
rivsites{ridx} = wtmp(isclose);

% Save table

RivList.sites = rivsites;
RivList = struct2table(RivList);

%% Read in full timeseries for all sites in RivList

fprintf('Reading river discharge and temperature data...\n');

if Opt.fetchdata
    riverdatadir = downloadusgsriverdata(RivList, 'format', 'json');
else
    rivdirs = dir('rivertimeseries*');
    riverdatadir = rivdirs(end).name; % Most recent one
end

Tmp = load(fullfile(riverdatadir, 'rivts.mat'), 'D');
D = Tmp.D;
    
%% Combine the USGS dataset with the Former Soviet Union datset
 
fprintf('Adding Former Soviet Union data to lists...\n');

% List of rivers of interest and their mouth coordinates

rusrivers = {...
    'Kamchatka' [56 12 28.17] [162 29  7.83]
    'Avacha'    [53  2  7.74] [158 30 42.62]
    'Anadyr'    [64 45  6.79] [176 25 36.19]};
RusRiv = cell2table(rusrivers, 'variableNames', {'name', 'lat', 'lon'});
RusRiv.lat = dms2degrees(RusRiv.lat);
RusRiv.lon = dms2degrees(RusRiv.lon);

nrus = height(RusRiv);

% Read .inv data (info on the stations)
% variable	format	width	1stcol	lastcol 

fmt = {...
'gauge no.'			'i5'	 5	  1	  5  
'blank'				'1x'	 1	  6	  6
'FSU gauge ID'		'i5'	 5	  7	 11
'blank'				'1x'	 1	 12	 12
'river name'		'a30'	30	 13	 42
'gauge name'		'a35'	35	 43	 77
'country code'		'a2'	 2	 78	 79
'blank'				'1x'	 1	 80	 80
'continent code'	'a1'	 1	 81	 81
'latitude'			'f8.3'	 8	 82  89
'longitude'			'f9.3'	 9	 90	 98
'drainage area'		'f11.1'	11	 99	109
'elevation'			'i5'	 5	110	114};

txt = fileread(Opt.fsuinvfile);
txt = regexp(txt, '\n', 'split');
txt = strvcat(txt{:});

idx = find(~strcmp(fmt(:,1), 'blank'));

isnum = regexpfound(fmt(:,2), '^[if]');

data = cell(1, length(idx));
for ii = 1:length(idx)
    data{ii} = cellstr(txt(3:end, fmt{idx(ii),4}:fmt{idx(ii),5}));
    if isnum(idx(ii))
        data{ii} = cellfun(@str2double, data{ii});
    end
end

vname = regexprep(fmt(idx,1), '[^A-Za-z_]', '');
RusInv = table(data{:}, 'VariableNames', vname);

% Read data

rusdata = load(Opt.fsudatafile); % m^3/s

isin = ismember(RusInv.rivername, RusRiv.name);

tf = ismember(rusdata(:,1), RusInv.gaugeno(isin));
rusdata = rusdata(tf,:);
rusdata(rusdata == -9) = NaN;

% Reformat to match datetime/value tables in D, and convert from m^3/s to
% ft^3/s to match USGS data

[gnum, data] = aggregate(rusdata(:,1), rusdata(:,2:end));
gdata = cell(nrus,1);

rname = vlookup(RusInv, gnum, 'rivername', 'gaugeno');
gname = vlookup(RusInv, gnum, 'gaugename', 'gaugeno');
Fsu.lat = vlookup(RusInv, gnum, 'latitude', 'gaugeno');
Fsu.lon = vlookup(RusInv, gnum, 'longitude', 'gaugeno');
Fsu.name = cellfun(@(a,b) sprintf('%s River: %s',a,b), rname, gname, 'uni', 0);
Fsu.code = cellstr(num2str(gnum, 'FSU-%d'));

m2ft = 100/2.54/12;

nusgs = height(D);
for ig = 1:length(gnum)
    
    yr = repmat(data{ig}(:,1), 1, 12);
    mn = repmat(1:12, size(yr,1), 1);
    dy = ones(size(yr)) * 15;
    
    t = datetime(yr(:), mn(:), dy(:));
    x = data{ig}(:,2:end)*(m2ft^3);
    
    isn = isnan(x(:)) | t < datetime(1950,1,1);
    t = t(~isn);
    x = x(~isn);
    
    [~,isrt] = sort(t);

    Fsu.data{ig,1} = table(t(isrt), x(isrt), 'variablenames', {'time','value'});

end

D = cat(1, D, struct2table(Fsu));

% Add the FSU rivers to the river list

[rivname, russites] = aggregate(rname, Fsu.code);

ruslat = vlookup(RusRiv, rivname, 'lat', 'name');
ruslon = vlookup(RusRiv, rivname, 'lon', 'name');

rivname = cellfun(@(x) sprintf('%s R', x), rivname, 'uni', 0);

RivNew = table(rivname, ruslat, ruslon, russites, ...
    'variablenames', {'name', 'lat', 'lon', 'sites'});

RivList = [RivList; RivNew];

%% Choose which stations to use for each river
% 
% For discharge, which rivers have multiple stations?  Which is closest to
% the mouth and most complete?  Is there any overlap in the data?  This
% step involves simply looking at each river, and taking notes that will be
% used in the next cell to build timeseries.

% Choose data for stations with multiple sites available

if Opt.choosedata

    Rworld = shapeprjread(Opt.nerivfile); % for Russian
    
    h.fig = figure;
    h.ax = axes('position', [0.05 0.05 0.95 0.9]);
    h.ax = subgrid(h.ax, 1, [0.7 0.3]); 
    h.ax = cat(2, h.ax{:});
    axes(h.ax(2));
    
    usamap(minmax(D.lat), minmax(wrapTo360(D.lon)));
    borders('russia', 'color', rgb('brown'));
    borders('alaska', 'color', rgb('brown'));
    
    is1 = [R.HEIRARCHY] == 1;
    is2 = [R.HEIRARCHY] == 2;

    plotm([R(is2).Lat], [R(is2).Lon], 'color', rgb('navy blue'));
    plotm([R(is1).Lat], [R(is1).Lon], 'color', rgb('navy blue'), 'linewidth', 2);
    plotm([Rworld.Lat], [Rworld.Lon], 'color', rgb('navy blue'), 'linewidth', 2);
    
    plotm(boxlat2, boxlon2, 'color', rgb('orange'));

    for ii = 1:nriv
        if length(RivList.sites{ii}) > 1
            
            title(h.ax(2), RivList.name{ii});
            
            ddata = vlookup(D, RivList.sites{ii}, 'data', 'code');
            isemp = cellfun('isempty', ddata); % station w/ T but not D
            sitetmp = RivList.sites{ii}(~isemp);
            
            npt = length(sitetmp);
            
            lt = [vlookup(D, sitetmp, 'lat', 'code') nan(npt,1)]';
            ln = [vlookup(D, sitetmp, 'lon', 'code') nan(npt,1)]';
            ddata = ddata(~isemp);
            
            xplt = cellfun(@(x) x.time, ddata, 'uni', 0);
            yplt = cellfun(@(x) x.value, ddata, 'uni', 0);
            
            for iy = 1:length(yplt)
                yplt{iy}(yplt{iy} == -999999) = NaN;
            end
            xy = [xplt yplt]';
            
            axes(h.ax(2));
            h.m = plotm(RivList.lat(ii), RivList.lon(ii), 'r*');
            h.s = plotm(lt, ln, 'x');
            set([h.m; h.s], 'markersize', 10, 'linewidth', 2);

            h.ts = plot(h.ax(1), xy{:});
            h.leg = legend(h.ts, sitetmp);

            pause;
            delete(h.m);
            delete(h.s);
            delete(h.ts);

        end
    end
    close(h.fig);
end

%% Build timeseries

fprintf('Building timeseries...\n');

% Put all the timeseries on the same time grid.  Not interpolating, just
% plugging in values where I have them.

tall = cellfun(@(x) x.time, D.data, 'uni', 0);
tall = unique(cat(1, tall{:}));
nt = length(tall);

discharge = nan(nt, size(D,1));
for ii = 1:size(D,1)
    [tf,loc] = ismember(D.data{ii}.time, tall);
    discharge(loc,ii) = D.data{ii}.value;
end
discharge(discharge == -999999) = NaN;

nriv = height(RivList);

%------------------------------
% Step 1: Fill in data for 
% rivers with only one station 
% available
%------------------------------

Dts.filled = nan(nt,nriv);  % filled-in timeseries data
Dts.filltype = zeros([size(Dts.filled) 2]); % track which process was used to fill, and which site

for ii = 1:nriv
    if length(RivList.sites{ii}) == 1
        tf = ismember(RivList.sites{ii}, D.code);
        if tf
            [Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, RivList.name{ii}, RivList.sites{ii}, RivList.name, D.code, Dts.filltype);
        end
    end
end

%------------------------------
% Step 2: For remaining rivers,
% apply data as noted in 
% previous cell
%------------------------------

% Note: switching from 'engine' to 'fit' in the slmfitrivers calls will
% bring up plots of the specific fits of overlapping data

% SLM prescription used as a base for all of them

prescrip = {'increasing', 'on', 'extrapolation', 'linear'};

% KOBUK R: Start with 15744500, then scale 15744000 to extend

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'KOBUK R', '15744500', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = slmfitrivers(discharge, Dts.filled, 'KOBUK R', '15744500', '15744000', 'engine', ...
    RivList.name, D.code, Dts.filltype, prescrip{:}, 'rightminslope', 1);

% STEWART R: Use 15625900

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'STEWART R', '15625900', RivList.name, D.code, Dts.filltype);

% YUKON R: Start with 15565447.  Fill in with 15564800 for early period
% (tested scaling, but the overlapping period doesn't capture high values,
% so I get better-matching variability just going with the upstream
% station).  For remaining gaps, scale 15453500.  Checked 15356000 as a
% scaling candidate, but that one's too far upstream... couldn't find a
% prescription that didn't lead to an mean/variablity shift in the
% resulting timeseries. 

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'YUKON R', '15565447', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'YUKON R', '15564800', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = slmfitrivers(discharge, Dts.filled, 'YUKON R', '15565447', '15453500', 'engine', ...
    RivList.name, D.code, Dts.filltype, prescrip{:});

% KUSKOKWIM R: Use 15304000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'KUSKOKWIM R', '15304000', RivList.name, D.code, Dts.filltype);

% WOOD R: Use 15303000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'WOOD R', '15303000', RivList.name, D.code, Dts.filltype);

% TERROR R: Use 15295700

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'TERROR R', '15295700', RivList.name, D.code, Dts.filltype);

% SUSITNA R: Use 15294350, scale 15292000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'SUSITNA R', '15294350', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = slmfitrivers(discharge, Dts.filled, 'SUSITNA R', '15294350', '15292000', 'engine', ...
    RivList.name, D.code, Dts.filltype, prescrip{:});

% SHIP C: Use 15276000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'SHIP C', '15276000', RivList.name, D.code, Dts.filltype);

% RESURRECTION C: Use both 15267900 and 15268000 (no overlap)

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'RESURRECTION C', '15267900', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'RESURRECTION C', '15268000', RivList.name, D.code, Dts.filltype);

% KENAI R: Use 15266300, scale 15258000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'KENAI R', '15258000', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = slmfitrivers(discharge, Dts.filled, 'KENAI R', '15266300', '15258000', 'engine', ...
    RivList.name, D.code, Dts.filltype, prescrip{:});

% ANCHOR R: Use both 15239900 and 15240000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'ANCHOR R', '15239900', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'ANCHOR R', '15240000', RivList.name, D.code, Dts.filltype);

% DUCK R: Use 15224000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'DUCK R', '15224000', RivList.name, D.code, Dts.filltype);

% LOWE R: Use 15226620, 15226600, and 15226500

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'LOWE R', '15226620', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'LOWE R', '15226600', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'LOWE R', '15226500', RivList.name, D.code, Dts.filltype);

% COPPER R: Use 15214000, scale 15212000

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'COPPER R', '15214000', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = slmfitrivers(discharge, Dts.filled, 'COPPER R', '15214000', '15212000', 'engine', ...
    RivList.name, D.code, Dts.filltype, prescrip{:});

% BATTLE C: Use 15238985 and 15238986

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'BATTLE C', '15238985', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled, 'BATTLE C', '15238986', RivList.name, D.code, Dts.filltype);

% LOWELL C: Use 1523849020 and 15238500

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled,  'LOWELL C', '1523849020', RivList.name, D.code, Dts.filltype);
[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled,  'LOWELL C', '15238500', RivList.name, D.code, Dts.filltype);

% CHESTER C: Use 15275100, scale 15275000 (this scaling needed a
% point-specific wieghting to throw out one overly-influential point, so
% doing it outside the function)  

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled,  'CHESTER C', '15275100', RivList.name, D.code, Dts.filltype);

[~,loc] = ismember({'15275100', '15275000'}, D.code);
[~,idx] = ismember('CHESTER C', RivList.name);
y = discharge(:, loc(1));
x = discharge(:, loc(2));
hasboth = ~isnan(x) & ~isnan(y);

[~,imax] = max(x(hasboth));
wght = ones(size(x(hasboth)));
wght(imax) = 0.1;

slm = slmengine(x(hasboth), y(hasboth), prescrip{:}, 'rightmaxslope', 1.5, 'weight', wght);
fill = isnan(Dts.filled(:,idx)) & ~isnan(x);
Dts.filled(fill,idx) = slmeval(x(fill), slm);
Dts.filltype(fill,idx,1) = 2;
Dts.filltype(fill,idx,2) = loc(2);

% Anadyr R: Use FSU-95051

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled,  'Anadyr R', 'FSU-95051', RivList.name, D.code, Dts.filltype);

% Kamchatka R: Use FSU-90997

[Dts.filled, Dts.filltype] = fillriver(discharge, Dts.filled,  'Kamchatka R', 'FSU-90997', RivList.name, D.code, Dts.filltype);

% ** Special handling **

% Note: BUSKIN R, HUMBOLT C, and STEWART R all have too little data to
% really create any sort of seasonal cycle.  We won't use these.

toolittle = ismember(RivList.name, {'BUSKIN R', 'HUMBOLT C', 'STEWART R'});

% Eklutna Creek: Upper Dam was built after this site goes inactive
% (decrease in flow is apparent in last few years of data).  Fill after
% this with zeros to prevent insertion of climatology in next step.

isek = strcmp(RivList.name, 'EKLUTNA C');
isn = isnan(Dts.filled(:,isek));
Dts.filled(isn,isek) = 0;
Dts.filltype(isn,isek) = 1;

%------------------------------
% Step 3: Build monthly 
% averages, interannual and 
% climatological
%------------------------------

[monthyear, Dts.monthly] = aggregate([year(tall) month(tall)], Dts.filled, @(x) nanmean(x,1));
Dts.monthly = cat(1, Dts.monthly{:});

tmonth = datetime([monthyear ones(size(monthyear,1),1)*15]);

[monthclima, Dts.clima] = aggregate(month(tall), Dts.filled, @(x) nanmean(x,1));
Dts.clima = cat(1, Dts.clima{:});

Dts.filltypemonth = zeros(size(Dts.monthly));
Dts.fillsitemonth = zeros(size(Dts.monthly));

for ir = 1:nriv
    [mnyr, tmp] = aggregate([year(tall) month(tall)], permute(Dts.filltype(:,ir,:), [1 3 2]));
    for itmp = 1:size(tmp,1)
        [tmp2, ntmp] = aggregate(tmp{itmp}, tmp{itmp}, @(x) size(x,1));
        [tmp2,ia] = setdiff(tmp2, [0 0], 'rows');
        if ~isempty(tmp2)
            ntmp = cell2mat(ntmp(ia));
            if size(tmp2,1) > 1
                [~,imax] = max(ntmp);
                Dts.filltypemonth(itmp,ir) = tmp2(imax,1);
                Dts.fillsitemonth(itmp,ir) = tmp2(imax,2);
            else
                Dts.filltypemonth(itmp,ir) = tmp2(1);
                Dts.fillsitemonth(itmp,ir) = tmp2(2);
            end
        end
    end
end

% [~, filltype] = aggregate([year(tall) month(tall)], Dts.filltype(:,:,1));
% [~, fillsite] = aggregate([year(tall) month(tall)], Dts.filltype(:,:,2));
% Dts.filltypemonth = zeros(size(Dts.monthly));
% Dts.fillsitemonth = cell(size(Dts.monthly));
% 
% return
% for ii = 1:length(filltype)
%     for ir = 1:size(filltype{ii},2)
%         unq = unique(filltype{ii}(:,ir));
%         if length(unq) == 1
%             Dts.filltypemonth(ii,ir) = unq;
%         else
%             tmp = filltype{ii}(:,ir);
%             tmp = tmp(tmp~=0);
%             unq = unique(tmp);
%             if length(unq) == 1
%                 Dts.filltypemonth(ii,ir) = unq;
%             else
%                 n = hist(tmp, unq);
%                 [~,imax] = max(n);
%                 Dts.filltypemonth(ii,ir) = unq(imax);
%             end
%         end
%     end
% end

%------------------------------
% Step 4: Fill interannual with
% climatology where missing
%------------------------------

[~,mloc] = ismember(monthyear(:,2), monthclima);
for ir = 1:nriv
    needfill = isnan(Dts.monthly(:,ir));
    Dts.monthly(needfill,ir) = Dts.clima(mloc(needfill),ir);
    Dts.filltypemonth(needfill,ir) = 3;
end

%------------------------------
% Step 5: Inpaint
%------------------------------

% for ir = 1:nriv
%     Dts.monthly(:,ir) = inpaint_nans(Dts.monthly(:,ir), 4);
% end

% Convert to kg/s, and eliminate any below-zero artifacts

cvt = 1000/(m2ft^3); % ft^3 -> kg, assuming fresh water
Dts.monthly = max(Dts.monthly, 0) * cvt; % kg/s

%% Distribute freshwater input across grid cells

fprintf('Distributing river input across grid cells...\n');

ell = referenceEllipsoid('earth');

[nxi, neta] = size(Grd.lat_rho);
nriv = height(RivList);

% Calculateriver influence on each grid cell based on e-folding scale

isin = inpolygon(wrapTo360(RivList.lon), RivList.lat, boxlon, boxlat);

frac = zeros(nxi, neta, nriv);

for ii = 1:nriv
    d = distance(Grd.lat_rho(:), Grd.lon_rho(:), RivList.lat(ii), RivList.lon(ii), ell)/1000;    
    if isin(ii)
    
        w = exp(-(1/Opt.efol) .* d);
        w(~Grd.mask_rho) = 0;
        wsum = sum(w);

        frac(:,:,ii) = reshape(w./wsum, nxi, neta);
    else % off grid
        frac(:,:,ii) = 0;
    end
        
end

% Grid cell area, in square meters

[ridx, cidx] = ndgrid(2:nxi-1, 2:neta-1);

rpsi = bsxfun(@plus, [-1 0 0 -1 -1], ridx(:));
cpsi = bsxfun(@plus, [-1 -1 0 0 -1], cidx(:));

rind = sub2ind([nxi, neta], ridx(:), cidx(:));
pind = sub2ind([nxi-1, neta-1], rpsi, cpsi);

plat = [Grd.lat_psi(pind) nan(numel(rind),1)]';
plon = [Grd.lon_psi(pind) nan(numel(rind),1)]';

grdarea = areaint(plat(:), plon(:), ell); % area, m^2

Grd.area = nan(nxi, neta);
Grd.area(2:end-1,2:end-1) = reshape(grdarea, nxi-2, neta-2);

% River flow will be added as extra precipitation to the ROMS model, in
% kg/m^2/s

meanflow = mean(Dts.monthly,1);
[srt, isrt] = sort(meanflow, 'descend');
isrt = setdiff(isrt, find(toolittle), 'stable'); % Remove a few

nt = length(tmonth);

precip = zeros(nxi, neta, nt);

for ii = isrt
    tmp = bsxfun(@times, frac(:,:,ii)./Grd.area, permute(Dts.monthly(:,ii), [2 3 1]));
    precip = precip + tmp;
end

% Convert time to netcdf-style reference date and time since arrays

basedate = datetime(1900,1,1,0,0,0);
tdur = tmonth - basedate;

%% Write to file

if Opt.writefile

    fprintf('Writing to file...\n');

    outfile = sprintf('%s.%d-%d.efol%dkm.updated%s.nc', outbase, year(minmax(tmonth)), Opt.efol, datestr(now, 'yymmdd'));

    Nc = struct('Name', '/', ... 
               'Format', 'classic', ...
               'Filename', outfile);

    % Set global attributes.

    Nc.Attributes(1).Name      = 'type';
    Nc.Attributes(1).Value     = 'ROMS forcing file';

    Nc.Attributes(2).Name      = 'title';
    Nc.Attributes(2).Value     = 'ROMS Runoff/precipitation for rivers';

    Nc.Attributes(3).Name      = 'grid_file';
    Nc.Attributes(3).Value     = grdfileshort;

    Nc.Attributes(4).Name      = 'base_date';
    Nc.Attributes(4).Value     = ['days since ' char(basedate)];

    Nc.Attributes(5).Name      = 'history';
    Nc.Attributes(5).Value     = sprintf('%s:  Created by K. Kearney with %s.m',     ...
                                        datestr(now),  mfilename);

    % Set file dimensions.

    Nc.Dimensions(1).Name      = 'xi_rho';
    Nc.Dimensions(1).Length    = nxi;
    Nc.Dimensions(1).Unlimited = false;

    Nc.Dimensions(2).Name      = 'eta_rho';
    Nc.Dimensions(2).Length    = neta;
    Nc.Dimensions(2).Unlimited = false;

    Nc.Dimensions(3).Name      = 'rain_time';
    Nc.Dimensions(3).Length    = nt;
    Nc.Dimensions(3).Unlimited = true;

    % Variables                                

    Nc.Variables(1) = roms_metadata('lat_rho',   true, 'nc_double', false);   
    Nc.Variables(2) = roms_metadata('lon_rho',   true, 'nc_double', false); 
    Nc.Variables(3) = roms_metadata('rain_time', true, 'nc_double', true); 
    Nc.Variables(4) = roms_metadata('rain',      true, 'nc_double', false); 

    Nc.Variables(4).Dimensions(1).Name = 'xi_rho';  % Change from lon
    Nc.Variables(4).Dimensions(2).Name = 'eta_rho'; % Change from lat

    % Fill

    Nc = check_metadata(Nc);

    % A few changes

    Nc.Variables(3).Attributes(2).Value = sprintf('days since %s', char(basedate));
    Nc.Variables(3).Attributes(1).Value = 'runoff time';
    Nc.Variables(4).Attributes(1).Value = 'runoff';
    Nc.Variables(4).Attributes(3).Value = 'runoff_time';
    Nc.Variables(4).Attributes(4).Value = 'xi_rho eta_rho runoff_time';

    % Write to file

    if exist(Nc.Filename, 'file')
        delete(Nc.Filename);
    end
    ncwriteschema(Nc.Filename, Nc);

    ncwrite(Nc.Filename, 'lat_rho', Grd.lat_rho);
    ncwrite(Nc.Filename, 'lon_rho', Grd.lon_rho);
    ncwrite(Nc.Filename, 'rain_time', days(tdur));
    ncwrite(Nc.Filename, 'rain', precip);

    % Rename rain to runoff (easier here than adding to roms_metadata)

    cmd = sprintf('ncrename -d rain_time,runoff_time -v rain,Runoff -v rain_time,runoff_time %s', Nc.Filename);
    system(cmd);

end







