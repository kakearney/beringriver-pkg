% This script covers the manual part of the river flow
% reconstruction... choosing which stations to use, and how to combine and
% interpolate when necessary. 
%
% The manual portion also includes adding in the 3 Russian rivers using the
% former Soviet Union dataset, and extending those timeseries via
% climatologies.
%
% The final variables of interest are:
%   Fuw:    table with metadata for Alaska rivers
%   RusRiv: table with metadata for Russia rivers
%   Ts:     structure with timeseries for all watersheds.  Columns
%           corresponds to the rows in Fuw and RusRiv.

%% Coastlines

if ~exist('cstflag', 'var')
    cstflag = false;
end

if cstflag

    % Lo-res coastlines (to get bounding box)

    [aklat, aklon] = borders('alaska');
    aklatlim = minmax(aklat);
    aklontmp = aklon;
    aklontmp(aklon > 0) = NaN;
    aklonlim = minmax(aklontmp);  % easier to ignore tail of Aleutians than combine

    % Hi-res coastlines (Note WDB = CIA World Data Bank, as opposed to WBD, watershed boundary dataset above)

    lev = 'f';
    cfile = fullfile(gshhsdir, 'GSHHS_shp', lev, sprintf('GSHHS_%s_L1.shp', lev));

    Cst = shapeprjread(cfile, 'BoundingBox', [aklonlim' aklatlim']);

end

%% Time
%
% Note: tested with both monthly and 5-day averaging.  My final analysis is
% a little hand-wavey... the results arer mostly the same, but I just
% prefer to go through the 5-day averaging first, then fill in via
% scaling/climatology, then doing a monthly average at the end.  This
% prevents short bits of a month (i.e. a day or two) from being used to
% represent an entire month when scaled/climatological data might be more
% appropriate.

yr = 1950:year(today);

% Short-term binning and averaging: 5-day averages (or 6-day, for leap day)

tbin = cell(length(yr),1);
dtbin = days(5);
for iyr = 1:length(yr)
    if leapyear(yr(iyr))
        tbin{iyr} = [datetime(yr(iyr),1,1):dtbin:datetime(yr(iyr),2,25) ...
               datetime(yr(iyr),3,2):dtbin:datetime(yr(iyr)+1,1,1)];     
    else
        tbin{iyr} = datetime(yr(iyr),1,1):dtbin:datetime(yr(iyr)+1,1,1);
    end
end
tbin = unique(cat(2, tbin{:}))';

tmid = mean([tbin(1:end-1) tbin(2:end)],2);

%% Set up arrays

dalaska = nan(length(tmid), height(Fuw)); % time x watershed array, river discharge, ft^3/s
datakeep = cell(height(Fuw), 1);

% return
% 
% [mn, yr] = ndgrid(1:12,1950:year(today));
% tbinmonthly = datetime(yr(:), mn(:), ones(numel(yr),1));
% tbinmonthly = tbinmonthly(tbinmonthly <= datetime('today'));
% 
% tbindaily = (datetime(1950,1,1):datetime('today'))';
% 
% tmid = mean([tbindaily(1:end-1) tbindaily(2:end)],2);


%% Watershed 1: 
%
% Lots of tributaries of the Copper R, which has two main stations:
% 15214000 is closer to the mouth and has a more recent record, but
% 15212000 includes most of the earlier period.  There are a few years of
% overlap, you we'll use that data to scale up the early period. Also add
% in a few small creeks that feed directly: 15216003 (MIDDLE ARM EYAK LK
% TR), 15216000 (POWER C), 15216008 (MURCHISON C), and 15215900 (GLACIER R
% TRIB).     

prescrip = {'increasing', 'on', 'extrapolation', 'linear'};
[ttmp, dtmp] = scaledata(Fuw.sitedata{1}, '15214000', '15212000', prescrip);
dtmp = max(dtmp, 0);

keep = ismember(Fuw.sitedata{1}.code, {'15216003','15216000','15216008','15215900'});

ttmp = [{ttmp}; Fuw.sitedata{1}.timed(keep)];
dtmp = [{dtmp}; Fuw.sitedata{1}.discharge(keep)];

dall = binandavg(tbin, ttmp, dtmp);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,1) = sum(dfill,2);

Fuw.sitedata{1}.used = keep | ismember(Fuw.sitedata{1}.code, {'15214000', '15212000'});


%% Watershed 2:
%
% There are two Lowell C stations, very close to each other and with no
% temporal overlap (one likely replaced the other).  Combine those two for
% Lowell C.  The others all represent individual streams; use all.

islowell = strcmp(Fuw.sitedata{2}.rivname, 'LOWELL C' );
ttmp = Fuw.sitedata{2}.timed(islowell);
dtmp = Fuw.sitedata{2}.discharge(islowell);
[ttmp, isrt] = sort(cat(1, ttmp{:}));
dtmp = cat(1, dtmp{:});
dtmp = dtmp(isrt);

ttmp = [{ttmp}; Fuw.sitedata{2}.timed(~islowell)];
dtmp = [{dtmp}; Fuw.sitedata{2}.discharge(~islowell)];

dall = binandavg(tbin, ttmp, dtmp);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,2) = sum(dfill, 2);

Fuw.sitedata{2}.used = true(height(Fuw.sitedata{2}),1);

%% Watershed 3
%
% Only two stations, and one only has one year worth of data.  Toss that
% one, use the more recent, longer running one.

keep = strcmp(Fuw.sitedata{3}.code, '15238648');

dall = binandavg(tbin, Fuw.sitedata{3}.timed(keep), Fuw.sitedata{3}.discharge(keep));

dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,3) = dfill;

Fuw.sitedata{3}.used = keep;

%% Watershed 4
%
% Lots of little rivers, several with multiple gauges.  The Solomon Gulch
% ones overlap for most time; 15226000 is most complete.  For Lowe R, the
% three don't overlap in time and are about the same magnitude, so combine.
% The Duck R also has overlapping time, use 15224000.  Finally, throw out
% Wolverine Glacier one; the Creek is further downstream and has more data.

islowe = strcmp(Fuw.sitedata{4}.rivname, 'LOWE R');
ttmp = Fuw.sitedata{4}.timed(islowe);
dtmp = Fuw.sitedata{4}.discharge(islowe);
[ttmp, isrt] = sort(cat(1, ttmp{:}));
dtmp = cat(1, dtmp{:});
dtmp = dtmp(isrt);

keep = ~islowe & ~ismember(Fuw.sitedata{4}.code, {'15225996', '15225997', '15225998', '15223900', '15236895'});

ttmp = [{ttmp}; Fuw.sitedata{4}.timed(keep)];
dtmp = [{dtmp}; Fuw.sitedata{4}.discharge(keep)];

dall = binandavg(tbin, Fuw.sitedata{4}.timed(keep), Fuw.sitedata{4}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,4) = sum(dfill,2);

Fuw.sitedata{4}.used = keep | islowe;

%% Watershed 5
%
% For the Terror R, there are two stations; 15295700 has more recent data,
% and we can fill in a few extra years by scaling the other.  Add 15296000
% (Uganik R) to that. 

prescrip = {'increasing', 'on', 'extrapolation', 'linear'};
[ttmp, dtmp] = scaledata(Fuw.sitedata{5}, '15295700', '15295600', prescrip);

keep = strcmp(Fuw.sitedata{5}.code, '15296000');
ttmp = [{ttmp}; Fuw.sitedata{5}.timed(keep)];
dtmp = [{dtmp}; Fuw.sitedata{5}.discharge(keep)];

dall = binandavg(tbin, ttmp, dtmp);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,5) = sum(dfill,2);

Fuw.sitedata{5}.used = ismember(Fuw.sitedata{5}.code, {'15295700', '15295600', '15296000'}); 

%% Watershed 6
%
% Throw out upstream tributaries of the Karluk R (15296550, 15296520,
% 15296500).  Use the rest. 

keep = ~ismember(Fuw.sitedata{6}.code, {'15296550', '15296520','15296500'});

dall = binandavg(tbin, Fuw.sitedata{6}.timed(keep), Fuw.sitedata{6}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,6) = sum(dfill,2);

Fuw.sitedata{6}.used = keep;

%% Watershed 7
%
% Two rivers, one station per river.  Not much data, but no decisions to be
% made here.

dall = binandavg(tbin, Fuw.sitedata{7}.timed, Fuw.sitedata{7}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,7) = sum(dfill,2);

Fuw.sitedata{7}.used = true(height(Fuw.sitedata{7}),1);

%% Watershed 8
%
% No decisions.  A few small rivers, one station per river. Use all.

dall = binandavg(tbin, Fuw.sitedata{8}.timed, Fuw.sitedata{8}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,8) = sum(dfill,2);

Fuw.sitedata{8}.used = true(height(Fuw.sitedata{8}),1);

%% Watershed 9
%
% Buskin R (15297437) has no data.  So left with just one river station, on
% Myrtle C (15297200)

keep = strcmp(Fuw.sitedata{9}.code, '15297200');

dall = binandavg(tbin, Fuw.sitedata{9}.timed(keep), Fuw.sitedata{9}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,9) = dfill;

Fuw.sitedata{9}.used = keep;

%% Watershed 10
%
% Two stations, both on Hidden Creek, both with the same one year worth of
% coverage.  Use the more downstream 15297110.

keep = strcmp(Fuw.sitedata{10}.code, '15297100');

dall = binandavg(tbin, Fuw.sitedata{10}.timed(keep), Fuw.sitedata{10}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,10) = dfill;

Fuw.sitedata{10}.used = keep;

%% Watershed 11
%
% Two stations. Alec R is tributary of Chignik R, so just use the Chignik
% one (15297585).

keep = strcmp(Fuw.sitedata{11}.code, '15297585');

dall = binandavg(tbin, Fuw.sitedata{11}.timed(keep), Fuw.sitedata{11}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,11) = dfill;

Fuw.sitedata{11}.used = keep;

%% Watershed 12
%
% Cook Inlet, land of all the rivers!  Including several with high
% streamflow (Susitna, Knik, Kenai, Matanuska, Chakachatna, Kasilof...)
%
% Clockwise from southwest... starting
% with some one station per stream ones: Paint R (15294900), Johnson R
% (15294700), Lake Fork Crecent R (15294640), Montana Bill C (15294585),
% Chakachatna R (15294500), Chuitna R (15294450).  Most downstream Susitna
% station is 15294350 (also has highest flow of all stations here); not
% very recent data but no great candidates for scaling.  Emptying into the
% Knik Arm are Lower Susitna (15290000), Cottonwood C (15286000), Wasilla C
% (15285000), Matanuska R (15284000), Knik R (15281000), Eklutna R
% (15280200), Peters C (15277410), Eagle R (15277100), Ship C (15276570),
% Chester C (15275100).  Continuing south, there Campbell C (15274600),
% Rabbit C (15273050). Turnagain Arm receives input from Resurrection C
% (15268000), Glacier C (15272550), Twentymile R (15272380), Portage C
% (15272280), Sixmile C (15271000). Continuing south... Bishop C
% (15267000), Kenai R (15266300), Kasilof R (15242000), Ninilchik R
% (15241600), Anchor R (15240000).  In Katchemak Bay, we have Fritz C
% (15239500), Bradley R (15239060... note, otpimizes for timeseries length
% and closeness to mouth), Battle C (15238986).  Continuing south, Barbara
% C (15238820), Seldovia R (15238795).  Didn't bother looking for scalable
% stuff, since most were available for small creeks that are negligible
% compared to the big rivers here.

keep = ismember(Fuw.sitedata{12}.code, {'15294900','15294700','15294640',...
    '15294585','15294500','15294450','15294350','15290000','15286000',...
    '15285000','15284000','15281000','15280200','15277410','15277100',...
    '15276570','15275100','15274600','15273050','15268000','15272550', ...
    '15272380','15272280','15271000','15267000','15266300','15242000', ...
    '15241600','15240000','15239500','15239060','15238820','15238795'});

dall = binandavg(tbin, Fuw.sitedata{12}.timed(keep), Fuw.sitedata{12}.discharge(keep));
dfill = max(0, inpaintn(fillclima(tbin(1:end-1), dall)));

dalaska(:,12) = sum(dfill,2);

Fuw.sitedata{12}.used = keep;

%% Watershed 13
%
% Two small creeks, not much data.  Use both.  There's not quite enough for
% the whole yearly climatology, so fill with inpainting (the inpaintn
% algorithm works better than inpaint_nan to maintain the "feel" of the
% dataset).

dall = binandavg(tbin, Fuw.sitedata{13}.timed, Fuw.sitedata{13}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dsum = sum(dfill,2);

dalaska(:,13) = inpaintn(dsum); 

Fuw.sitedata{13}.used = true(height(Fuw.sitedata{13}),1);

%% Watershed 14
%
% Only one station... use it

dall = binandavg(tbin, Fuw.sitedata{14}.timed, Fuw.sitedata{14}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,14) = dfill;

Fuw.sitedata{14}.used = true(height(Fuw.sitedata{14}),1);

%% Watershed 15
%
% Lots of little creeks.  This is an island HU, so no graph network...
% everything dumps into the ocean.  So use everything and sum.

dall = binandavg(tbin, Fuw.sitedata{15}.timed, Fuw.sitedata{15}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,15) = sum(dfill,2);

Fuw.sitedata{15}.used = true(height(Fuw.sitedata{15}),1);

%% Watershed 16
%
% Most stations are upstream of the Kvichak R (with most of their
% streamflow dumped into a lake just north of the Kvichak station).  So use
% only the Kvichak (15300500) and Eskimo C (15297900).

keep = ismember(Fuw.sitedata{16}.code, {'15297900', '15300500'});

dall = binandavg(tbin, Fuw.sitedata{16}.timed(keep), Fuw.sitedata{16}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,16) = sum(dfill,2);

Fuw.sitedata{16}.used = keep;

%% Watershed 17
%
% Three main non-tribs: Nushagak (15302500), Wood R (15303000), and Snake R
% (15303150).  The Nuyakuk (15302000) is a tributary of the Nushagak that
% has a pretty complete record; scale and use this one to extend the
% Nushagak record.

prescrip = {'increasing', 'on', 'extrapolation', 'linear'};
[ttmp, dtmp] = scaledata(Fuw.sitedata{17}, '15302500', '15302000', prescrip);

keep = ismember(Fuw.sitedata{17}.code, {'15303000', '15303150'});
ttmp = [{ttmp}; Fuw.sitedata{17}.timed(keep)];
dtmp = [{dtmp}; Fuw.sitedata{17}.discharge(keep)];

dall = binandavg(tbin, ttmp, dtmp);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,17) = sum(dfill,2);

Fuw.sitedata{17}.used = keep | ismember(Fuw.sitedata{17}.code, {'15302500', '15302000'});

%% Watershed 18
%
% Only one station with about a year worth of data.  Guess we'll use it.

dall = binandavg(tbin, Fuw.sitedata{18}.timed, Fuw.sitedata{18}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,18) = dfill;

Fuw.sitedata{18}.used = true(height(Fuw.sitedata{18}),1);


%% Watershed 19
%
% The Kuskokwim station near Crooked C (15304000) has a near-complete
% record. A few small tributaries are located downstream of this one
% (Browns C, Kisalarek R), but their streamflow is negligible compared to
% the Kuskokwim.

keep = strcmp(Fuw.sitedata{19}.code, '15304000');

dall = binandavg(tbin, Fuw.sitedata{19}.timed(keep), Fuw.sitedata{19}.discharge(keep));
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,19) = sum(dfill,2);

Fuw.sitedata{19}.used = keep;


%% Watershed 20
%
% Only the Unalakleet has a gauge, so use it.

dall = binandavg(tbin, Fuw.sitedata{20}.timed, Fuw.sitedata{20}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,20) = dfill;

Fuw.sitedata{20}.used = true(height(Fuw.sitedata{20}),1);


%% Watershed 21
%
% Two stations on Stewart R; throw out the upstream one (15625850).  The
% rest represent individual creeks.  Some datasets are intermittent, so
% inpaint where climatology is missing.

keep = ~ismember(Fuw.sitedata{21}.code, {'15625850'});

dall = binandavg(tbin, Fuw.sitedata{21}.timed(keep), Fuw.sitedata{21}.discharge(keep));
dfill = max(inpaintn(fillclima(tbin(1:end-1), dall)), 0);

dalaska(:,21) = sum(dfill,2);

Fuw.sitedata{21}.used = keep;


%% Watershed 22
%
% One gauge on El Dorado Creek.  Missing data during low-flow times, so
% inpaint there.

dall = binandavg(tbin, Fuw.sitedata{22}.timed, Fuw.sitedata{22}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,22) = max(inpaintn(dfill), 0);

Fuw.sitedata{22}.used = true(height(Fuw.sitedata{22}),1);

%% Watershed 23
%
% Two gauged streams... use them both.

dall = binandavg(tbin, Fuw.sitedata{23}.timed, Fuw.sitedata{23}.discharge);
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,23) = sum(dfill,2);

Fuw.sitedata{23}.used = true(height(Fuw.sitedata{23}),1);


%% Watershed 24
%
% Use the downstream Kobuk R (15744500), Noatak (15746000) and June C
% (15743000).  Humbolt C has too little data, and the rest are
% upstream/tributary Kobuk ones.


keep = ismember(Fuw.sitedata{24}.code, {'15744500', '15746000', '15743000'});

dall = binandavg(tbin, Fuw.sitedata{24}.timed(keep), Fuw.sitedata{24}.discharge(keep));
dfill = max(inpaintn(fillclima(tbin(1:end-1), dall)), 0);

dalaska(:,24) = sum(dfill,2);

Fuw.sitedata{24}.used = keep;

%% Watershed 25
%
% The Yukon!  15565447 is the most downstream station; everything else is
% upstream/a tributary. 15453500 is next closest to the mouth, and can be
% scaled to fill a gap in the 90s.  The most complete record available
% comes from the most upstream Yukon station (15356000); can be scaled to
% fill earlier years.

prescrip = {'increasing', 'on', 'extrapolation', 'linear'};
[ttmp1, dtmp1] = scaledata(Fuw.sitedata{25}, '15565447', '15453500', prescrip);
[ttmp2, dtmp2] = scaledata(Fuw.sitedata{25}, '15565447', '15356000', prescrip);
[ttmp3, i3] = setdiff(ttmp2, ttmp1);
[ttmp, isrt] = sort([ttmp1; ttmp3]);
dtmp = [dtmp1; dtmp2(i3)];
dtmp = dtmp(isrt);

dall = binandavg(tbin, {ttmp}, {dtmp});
dfill = fillclima(tbin(1:end-1), dall);

dalaska(:,25) = dfill;

Fuw.sitedata{25}.used = ismember(Fuw.sitedata{25}.code, {'15565447', '15453500', '15356000'});

%% Calculate monthly averages for all Alaska data

[yy,mm] = ndgrid(yr(1):yr(end), 1:12);
Ts.tbin = datetime(yy(:), mm(:), ones(numel(yy),1));
Ts.tbin = [sort(Ts.tbin); datetime(yr(end)+1,1,1)];

nt = length(Ts.tbin) - 1;

binidx = discretize(tmid, Ts.tbin);
[idxtmp, aktmp] = aggregate(binidx, dalaska, @(x) mean(x,1));

[~,loc] = ismember(idxtmp, 1:length(Ts.tbin)-1);
Ts.alaska = nan(nt, nwatershed);
Ts.alaska(loc,:) = cat(1, aktmp{:});

%% Build timeseries for Russian rivers

% List of rivers of interest, their mouth coordinates, and the gauge number
% associated with the most downstream station for each

rusrivers = {...
    'Kamchatka' [56 12 28.17] [162 29  7.83]    90997
    'Avacha'    [53  2  7.74] [158 30 42.62]    90926
    'Anadyr'    [64 45  6.79] [176 25 36.19]    95051};
RusRiv = cell2table(rusrivers, 'variableNames', {'name', 'lat', 'lon', 'gaugeno'});
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

txt = fileread('/Volumes/Storage/UCAR_553p2_russiaRivers/fsu.inv');
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

rusdata = load('/Volumes/Storage/UCAR_553p2_russiaRivers/fsu2.txt'); % m^3/s

tf = ismember(rusdata(:,1), RusRiv.gaugeno);
rusdata = rusdata(tf,:);
rusdata(rusdata == -9) = NaN;

% Reformat into complete monthly timeseries, filling with climatology where
% necessary, and convert from m^3/s to ft^3/s to match USGS data

[gnum, data] = aggregate(rusdata(:,1), rusdata(:,2:end));

[~,loc] = ismember(RusRiv.gaugeno, gnum);
data = data(loc);

Ts.russia = nan(nt, 3); % Discharge array 
m2ft = 100/2.54/12;

RusRiv.t = cell(3,1);
RusRiv.x = cell(3,1);
for ii = 1:3
    yr = repmat(data{ii}(:,1), 1, 12);
    mn = repmat(1:12, size(yr,1), 1);
    dy = ones(size(yr)) * 15;
    
    RusRiv.t{ii} = datetime(yr(:), mn(:), dy(:));
    RusRiv.x{ii} = data{ii}(:,2:end)*(m2ft^3);
    
    idx = discretize(RusRiv.t{ii}, Ts.tbin);
    isn = isnan(idx);
    
    Ts.russia(idx(~isn),ii) = RusRiv.x{ii}(~isn);
    
end
Ts.russia = fillclima(Ts.tbin(1:end-1), Ts.russia);

%% Plot

h = plotgrid('function', {@(y) stairs(Ts.tbin,[y; NaN]/1000), num2cell([Ts.alaska Ts.russia],1)'}, ...
    'sp', 0, 'outputs', {'ln'}, 'ml', 0.2', 'mr', 0.02);
h.ln = cat(1, h.ln{:});
set(h.fig, 'color', 'w');
set(h.ax, 'box', 'off', 'ylim', [0 100], 'xlim', minmax(tmid), 'clipping', 'off');
set(h.ax(1:end-1), 'xcolor', 'none');
set(h.ax(2:2:end), 'yaxisloc', 'right');
h.lbl = multitextloc(h.ax, [Fuw.name; RusRiv.name], 'northwestoutside');
set(h.lbl, 'fontsize', 8);

cmap = cptcmap('Set2_08');
cmap = repmat(cmap, 4, 1);
cmap = num2cell(cmap(1:numel(h.ax),:),2);
set(h.ln, {'color'}, cmap);
set(h.ax, {'ycolor'}, cmap);
set(h.lbl, {'color'}, cmap);


%% Subfunctions


function x = binandavg(tbin, tcell, xcell)

    nbin = length(tbin) - 1;
    nx = length(xcell);
    x = nan(nbin,nx);
    for ii = 1:nx
        isn = xcell{ii} == -999999;
        tcell{ii} = tcell{ii}(~isn);
        xcell{ii} = xcell{ii}(~isn);
        tidx = discretize(tcell{ii}, tbin);
        [tidx, xavg] = aggregate(tidx, xcell{ii}, @mean);
        x(tidx,ii) = cat(1, xavg{:});
    end
end

function x = fillclima(t, x)
    dv = datevec(t);
    [mndy, data] = aggregate(dv(:,2:3), x, @(y) nanmean(y,1));
    data = cat(1, data{:});
    [~,loc] = ismember(dv(:,2:3), mndy, 'rows');
    data = data(loc,:);
    
    isn = isnan(x);
    x(isn) = data(isn);
   
end

function [ttmp, dtmp] = scaledata(F, codekeep, codescale, prescrip, flag)

    if nargin < 5
        flag = false;
    end
    
    isa = strcmp(F.code, codekeep);  % keep intact
    isb = strcmp(F.code, codescale); % to be scaled
    
    ta = F.timed{isa};
    tb = F.timed{isb};
    a = F.discharge{isa};
    b = F.discharge{isb};
    
    isn = a == -999999;
    a = a(~isn);
    ta = ta(~isn);
    
    isn = b == -999999;
    b = b(~isn);
    tb = tb(~isn);
    
    [tboth, ia, ib] = intersect(ta, tb);
    aoverlap = a(ia);
    boverlap = b(ib);
    
    if flag
        slm = slmfit(boverlap,aoverlap,prescrip{:});
        ttmp = {ta, tb};
        dtmp = {a, b};
        
    else
        slm = slmengine(boverlap,aoverlap,prescrip{:});

        [tbscale, ibscale] = setdiff(tb, ta);

        ttmp = [ta; tbscale];
        dtmp = [a; slmeval(b(ibscale), slm)];
        [ttmp, isrt] = sort(ttmp);
        dtmp = dtmp(isrt);
    end

end

