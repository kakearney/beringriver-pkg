function rivdatadir = downloadusgsriverdata(RivList, varargin)
%DOWNLOADUSGSRIVERDATA Read data from NWIS database
%
% riverdatadir = downloadusgsriverdata(RivList)
%
% This function downloads discharge and temperature data from the NWIS
% daily values database.  Once downloaded, the data is parsed from the
% individual ,xml files and saved to a .mat file as two tables.
%
% All data between 1950/01/01 and the present will be queried.
%
%
% Input variables:
%
%   RivList:        table of river names and their corresponding site
%                   numbers in the NWIS database
%
% Optional input variables, passed as parameter/value pairs:
%
%   format:         format to use when querying.  Can be 'waterml,1.1',
%                   'json', or 'waterml,2.0'.  
%                   Note: Currently (12/15/2016), the JSON format reflect
%                   waterML 1.1.  Looks like that may change in the future
%                   to mimic 2.0 instead.  If parsing fails, check this.
%                   Note 2: waterml,1.1 is now deprecated on the USGS
%                   website, so it can't be used to fetch data anymore.
%                   But it can still be used to rerun the parsing of data
%                   with pre-existing files.
%
%   folder:         folder where new queried data will be saved.  If folder
%                   already includes a file with a name matching that of
%                   one of the rivers to be queried, that river will be
%                   skipped in the fetch data portion of this function.
%
% Output variables:
%
%   riverdatadir:   folder where all downloaded .xml files, as well as the
%                   the rivts.mat file (holding all data) is saved.  Unless
%                   of folder is passed as input, will follow the pattern
%                   of rivertimeseriesYYYYMMDD, where YYYYMMDD is today's
%                   date.

% Copyright 2015-2016 Kelly Kearney

% Parse input

tnow = datestr(now, 'yyyy-mm-dd');

p = inputParser;
p.addParameter('format', 'waterml,1.1', @(x) validateattributes(x, {'char'}, {}));
p.addParameter('folder', sprintf('rivertimeseries%s', datestr(now, 'yyyymmdd')), @(x) validateattributes(x, {'char'}, {}));
p.parse(varargin{:});

Opt = p.Results;

rivdatadir = Opt.folder;

% Read site data, river by river (originally attempted all sites at once,
% but that was too big and kept timing out)

nriv = size(RivList,1);

if ~exist(rivdatadir, 'dir')
    mkdir(rivdatadir);
end
tsfile = cell(nriv,1);
for ir = 1:nriv
    
    switch Opt.format
        case {'waterml,1.1', 'waterml,2.0'}
            tsfile{ir} = fullfile(rivdatadir, sprintf('rivts_%s.xml', regexprep(RivList.name{ir}, '\s', '_')));
        case 'json'
            tsfile{ir} = fullfile(rivdatadir, sprintf('rivts_%s.json', regexprep(RivList.name{ir}, '\s', '_')));
    end
           
    sitestr = sprintf('%s,', RivList.sites{ir}{:});
    sitestr = sitestr(1:end-1);
    
    urlcmd = sprintf('https://waterservices.usgs.gov/nwis/dv/?format=%s&sites=%s&startDT=1950-01-01&endDT=%s&statCd=00003&parameterCd=00060,00010', Opt.format, sitestr, tnow);

    if ~exist(tsfile{ir}, 'file')
        fprintf('%s\n', RivList.name{ir});
        [f, status] = urlwrite(urlcmd, tsfile{ir});
    end
end


for ir = nriv:-1:1
    switch Opt.format
        case 'waterml,1.1'
            [D(ir), T(ir)] = parsewaterml1p1(tsfile{ir});
        case 'json'
            [D(ir), T(ir)] = parsejson(tsfile{ir});
    end
end

D = arrayfun(@struct2table, D, 'uni', 0);
D = cat(1, D{:});

T = arrayfun(@struct2table, T, 'uni', 0);
T = cat(1, T{:});

save(fullfile(rivdatadir, 'rivts'), 'D', 'T');

%--------------------------------
% Read data from waterML 1.0 file
%--------------------------------

function [D, T] = parsewaterml1p1(file)

try
    dom = xmlread(file);
catch
    tmp = fileread(file);
    error('Unable to parse .xml file; contents are as follows:\n\n %s', tmp);
end

% Read variable codes

tmp = dom.getElementsByTagName('ns1:variableCode');
nts = tmp.getLength;
vcode = cell(nts,1);
for ii = 1:nts
    vcode{ii} = char(tmp.item(ii-1).item(0).getData);
end

% Read discharge values

idx = find(strcmp(vcode, '00060')) - 1;
n = length(idx);

[D.name, D.code, D.data] = deal(cell(n,1));
[D.lat, D.lon] = deal(nan(n,1));

for ii = 1:n
    sinfo = tmp.item(idx(ii)).getParentNode.getParentNode;
    D.name{ii} = char(sinfo.getElementsByTagName('ns1:siteName').item(0).item(0).getData);
    D.code{ii} = char(sinfo.getElementsByTagName('ns1:siteCode').item(0).item(0).getData);
    D.lat(ii) = str2double(char(sinfo.getElementsByTagName('ns1:latitude').item(0).item(0).getData));
    D.lon(ii) = str2double(char(sinfo.getElementsByTagName('ns1:longitude').item(0).item(0).getData));

    nodata = str2double(char(sinfo.getElementsByTagName('ns1:noDataValue').item(0).item(0).getData));

    val = sinfo.getElementsByTagName('ns1:value');
    nv = val.getLength;
    [ttmp,dtmp] = deal(cell(nv,1));
    for iv = 1:nv
        ttmp{iv} = char(val.item(iv-1).getAttribute('dateTime'));
        dtmp{iv} = char(val.item(iv-1).item(0).getData);
    end
    ttmp = datetime(ttmp, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'format', 'dd-MM-uuuu');
    dtmp = str2num(strvcat(dtmp{:}));
    dtmp(dtmp == nodata) = NaN;
    D.data{ii} = table(ttmp, dtmp, 'VariableNames', {'time', 'value'});
end

% Read temperature values

idx = find(strcmp(vcode, '00010')) - 1;
n = length(idx);

[T.name, T.code, T.data] = deal(cell(n,1));
[T.lat, T.lon] = deal(nan(n,1));

for ii = 1:n
    sinfo = tmp.item(idx(ii)).getParentNode.getParentNode;
    T.name{ii} = char(sinfo.getElementsByTagName('ns1:siteName').item(0).item(0).getData);
    T.code{ii} = char(sinfo.getElementsByTagName('ns1:siteCode').item(0).item(0).getData);
    T.lat(ii) = str2double(char(sinfo.getElementsByTagName('ns1:latitude').item(0).item(0).getData));
    T.lon(ii) = str2double(char(sinfo.getElementsByTagName('ns1:longitude').item(0).item(0).getData));

    nodata = str2double(char(sinfo.getElementsByTagName('ns1:noDataValue').item(0).item(0).getData));

    val = sinfo.getElementsByTagName('ns1:value');
    nv = val.getLength;
    [ttmp,dtmp] = deal(cell(nv,1));
    for iv = 1:nv
        ttmp{iv} = char(val.item(iv-1).getAttribute('dateTime'));
        dtmp{iv} = char(val.item(iv-1).item(0).getData);
    end
    ttmp = datetime(ttmp, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'format', 'dd-MM-uuuu');
    dtmp = str2num(strvcat(dtmp{:}));
    dtmp(dtmp == nodata) = NaN;
    T.data{ii} = table(ttmp, dtmp, 'VariableNames', {'time', 'value'});
end

%--------------------------------
% Read data from JSON file
%--------------------------------

function [D, T] = parsejson(file)

try
    data = loadjson(file, 'SimplifyCell', 1);

    n = length(data.value.timeSeries);

    varcode = arrayfun(@(X) X.variable.variableCode.value,                   data.value.timeSeries, 'uni', 0)';
    DT.name = arrayfun(@(X) X.sourceInfo.siteName,                           data.value.timeSeries, 'uni', 0)';
    DT.code = arrayfun(@(X) X.sourceInfo.siteCode.value,                     data.value.timeSeries, 'uni', 0)';
    DT.lat  = arrayfun(@(X) X.sourceInfo.geoLocation.geogLocation.latitude,  data.value.timeSeries)';
    DT.lon  = arrayfun(@(X) X.sourceInfo.geoLocation.geogLocation.longitude, data.value.timeSeries)';

    DT.data = cell(size(varcode));
    for ii = 1:n
        val =  cellfun(@str2double, {data.value.timeSeries(ii).values.value.value}');
        ttmp = {data.value.timeSeries(ii).values.value.dateTime}';
        ttmp = datetime(ttmp, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'format', 'dd-MM-uuuu');

        DT.data{ii} = table(ttmp, val, 'VariableNames', {'time', 'value'});
    end

    fld = fieldnames(DT);
    isd = strcmp(varcode, '00060');

    for ii = 1:length(fld)
        D.(fld{ii}) = DT.(fld{ii})(isd);
        T.(fld{ii}) = DT.(fld{ii})(~isd);
    end
catch ME
    error('Failed to parse... did the JSON format update to 2.0?\n%s', ME.message);
end



