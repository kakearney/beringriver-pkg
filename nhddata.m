% This script downloads data for rivers along the coast of Alaska from the
% National Hydrography Database. It then extracts all the flowline data
% from the downloaded shapefiles and saves them to a .mat file.
%
% cartExport.csv file was created via the NHD interactive viewer, exported
% to csv.  It lists all the hydrographic units in my area of interest
% (coastal Alaska).

%% Extract ftp addresses for files

nhdfile = '/Volumes/Storage/NationalHydrographyDataset/cartExport.csv';
NhdInfo = readtext(nhdfile, ',', '', '"');
NhdInfo = cell2table(NhdInfo(2:end,:), 'VariableNames', NhdInfo(1,:));

[pth,fl,ex] = fileparts(NhdInfo.URL{1});
pthparts = regexp(pth, '/+', 'split');

serverfolder = fullfile(pthparts{3:end});

nhddir = '/Volumes/Storage/NationalHydrographyDataset/';


%% Fetch all flowline files via ftp

% Remove any failed download attempts

Tmp = dirfull(fullfile(nhddir, '*.zip'));
issmall = [Tmp.bytes] < (1024^2); % less than 1 MB
cellfun(@delete, {Tmp(issmall).name});

% Fetch data (note: this doesn't work, due to setup of their ftp server...
% had to do it via web manually).

f = ftp(pthparts{2});
cd(f, serverfolder);

for ii = 1:height(NhdInfo)
    [~,fl,ex] = fileparts(NhdInfo.URL{ii});
    ftpfile = [fl ex];
    
    fprintf('  Downloading %s\n', ftpfile);
    if ~exist(fullfile(nhddir, ftpfile), 'file')
        try
            mget(f, ftpfile, nhddir);
        catch ME
            fprintf('    Failed: %s\n', ME.message);
        end
    else
        fprintf('    Already done\n');
    end
end
close(f);

%% Extract rivers from NHD dataset

% NHD list of Alaskan features: extract streams

akfeaturefile = '/Volumes/Storage/NationalHydrographyDataset/AK_Features_20151001.txt';
fid = fopen(akfeaturefile);
header = fgetl(fid);
fclose(fid);
header = regexprep(header, '[^A-Z_|]', ''); % some weird non-ascii characters
header = regexp(header, '\|', 'split');
Nhd = readtable(akfeaturefile, 'filetype', 'text', 'delimiter', '|', ...
    'headerlines', 1, 'readvariablenames', false);
Nhd.Properties.VariableNames = header;
isriv = strcmp(Nhd.FEATURE_CLASS, 'Stream');
Nhd = Nhd(isriv,:);

%% Unzip NHD HU8 files and extract flow

Files = dirfull(fullfile(nhddir, '*HU8.zip'));
nfile = length(Files);

if exist('Shape', 'dir')
    rmdir('Shape', 's');
end
    
% Extract lat and lon of named rivers

[names, lat, lon] = deal(cell(nfile,1));

for ii = 1:length(Files)
    unzip(Files(ii).name);
    
    shpname = dirfull(fullfile('Shape', '*Flowline.shp'));

    if isempty(shpname) 
        warning('No flowline file found for %s', Files(ii).name);
    elseif length(shpname) > 1
        warning('Multiple possible flowline files for %s', Files(ii).name);
    else
        
        Data = shapeprjread(shpname.name);

        [gnisname, Data] = aggregate({Data.GNIS_NAME}', Data);
        isemp = cellfun('isempty', gnisname);

        names{ii} = gnisname(~isemp);
        lat{ii} = cellfun(@(X) [X.Lat], Data(~isemp), 'uni', 0);
        lon{ii} = cellfun(@(X) [X.Lon], Data(~isemp), 'uni', 0);
    end

    rmdir('Shape', 's');
    
end

isemp = cellfun('isempty', lat);
lat = lat(~isemp);
lon = lon(~isemp);

[ltpoly, lnpoly] = deal(cell(size(lat)));
for ii = 1:length(lat)

    lttmp = cat(2, lat{ii}{:});
    lntmp = cat(2, lon{ii}{:});
    lttmp = lttmp(~isnan(lttmp));
    lntmp = lntmp(~isnan(lntmp));
    k = convhull(lntmp, lttmp);

    ltpoly{ii} = lttmp(k);
    lnpoly{ii} = lntmp(k);

end

ltall = ltpoly{1};
lnall = lnpoly{1};
for ii = 2:length(lat)
    [ltall, lnall] = polybool('|', ltall, lnall, ltpoly{ii}, lnpoly{ii});
end
    
lat = cat(1, lat{:});
lon = cat(1, lon{:});
names = cat(1, names{:});

[lon, lat] = cellfun(@(a,b) joinsegments(a,b), lon, lat, 'uni', 0);
[lon, lat] = cellfun(@(x,y) polyjoin(x,y), lon, lat, 'uni', 0);

save nhddata lat lon names ltall lnall;


