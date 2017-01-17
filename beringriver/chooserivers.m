function chooserivers(matfile, S, R, Nhd, boxlat, boxlon, sname)
%CHOOSERIVERS Interactive plot to choose rivers
%
% chooserivers(matfile, S, R, Nhd, boxlat, boxlon, sname)
%
% This function plots a map of Alaska with lots of rivers (all big, and
% near-coast small ones), with an overlay of all the USGS stations and
% their names.  
%
% Left click to zoom in, double click to zoom out.  When you find a river
% that has a station on it and feeds into the ROMS domain (orange box),
% right click as close to the river mouth as you can.  An input box will
% prompt you for the name of the river you just clicked on.  When done,
% press enter.
%
% Input variables:
%
%   matfile:    name of file where river coordinates will be saved
%
%   S:          table array listing USGS station locations
%
%   R:          shapefile structure of major rivers
%
%   Nhd:        structure of National Hydrography Database rivers (see
%               nhddata.m)
%
%   boxlat:     latitude coordinates of box around model grid
%
%   boxlon:     longitude coordinates of box around model grid
%
%   sname:      river names, corresponding to rows in S

% Copyright 2015-2106 Kelly Kearney

figure('color', 'w');
axes('position', [0 0 1 1]);
usamap('Alaska');

is1 = [R.HEIRARCHY] == 1;
is2 = [R.HEIRARCHY] == 2;

[nhdlat, nhdlon] = polyjoin(Nhd.lat, Nhd.lon);
plotm(nhdlat, nhdlon, 'color', rgb('medium blue'));
plotm([R(is2).Lat], [R(is2).Lon], 'color', rgb('navy blue'));
plotm([R(is1).Lat], [R(is1).Lon], 'color', rgb('navy blue'), 'linewidth', 2);
borders('alaska', 'color', rgb('brown'));
plotm(boxlat, boxlon, 'color', rgb('orange'));

plotm(S.lat, S.lon, 'rx');
textm(S.lat, S.lon, sname);

mlat = nan(0,1);
mlon = nan(0,1);
mname = cell(0,1);
count = zeros;

while 1
    [x,y] = ginput2(1, 'KeepZoom');
    if ~isempty(x)
        count = count + 1;
        [mlat(count), mlon(count)] = minvtran(x,y);
        tmp = inputdlg('Name?');
        mname{count} = tmp{1};
    else
        break
    end
end

save(matfile,  'mlat',  'mlon', 'mname'); 
