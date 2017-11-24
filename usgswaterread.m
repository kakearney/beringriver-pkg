function Data = usgswaterread(file, format)
%USGSWATERREAD Read selected data from USGS Water Services data
%
% Data = usgswaterread(file)
% Data = usgswaterread(file, format)
%
% This functions parses data from files downloaded from the USGS Daily
% Values Site web service.  At the moment, it's tailored to extract only
% temperature and discharge data.
%
% Input variables:
%
%   file:   name of data file.
%
%   format: file format, can be either 'waterml,1.1' (WaterML, v1.1, the
%           older, now-deprecated XML format), 'waterml,2.0' (WaterML,
%           v2.0, the current XML format), or 'json' (JSON).  Default is
%           'waterml,1.1' if the input file has a .xml extension and 'json'
%           if the input file has a .json extension.
%
% Output variables:
%
%   Data

if nargin < 2
    [~,~,ex] = fileparts(file);
    if strcmpi(ex, '.json')
        format = 'json';
    elseif strcmpi(ex, '.xml')
        format = 'waterml,1.1';
    else
        error('Could not determine format from file extension');
    end
end
    

switch format
    case 'waterml,1.1'
        [D, T] = parsewaterml1p1(file);
    case 'json'
        [D, T] = parsejson(file);
    case 'waterml,2.0'
        error('watwerml,2.0 option not coded yet, sorry');
end

Data.discharge  = struct2table(D);
Data.temperature = struct2table(T);

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