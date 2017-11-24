% riversToRunoff
%
% This script translates the river runoff timeseries into spatially
% distributed freshwater flux.  The timeseries for the eastern side of the
% domain come from USGS data; each timeseries is associated with a single
% frontal hydrologic unit polygon (see riversfromwatersheds.m).  On the
% western side of the domain, we have three Russian rivers with monthly
% timeseries from an older former-Soviet dataset.

%% Run code for USGS data

watersheds;
riversfromwatersheds;

Opt.efol = 20;

%% Distribute freshwater input across grid cells

ell = referenceEllipsoid('earth');

[nxi, neta] = size(Grd.lat_rho);

% Calculate river influence on each grid cell based on e-folding scale

[lnmask, ltmask] = mask2poly(Grd.lon_psi, Grd.lat_psi, Grd.mask_rho(2:end-1,2:end-1));

nak = size(Ts.alaska,2);
nrus = size(Ts.russia,2);

frac = zeros(nxi, neta, nak+nrus);

for ii = 1:nak
    fprintf('Calculating distances: Alaska %d/%d\n', ii, nak);
    
    % Distance from each grid point to the polygon
    
    ltpoly = Wbd(Fuw.mouthidx(ii)).Lat;
    lnpoly = wrapTo360(Wbd(Fuw.mouthidx(ii)).Lon);
    
    xy = distance2curve([lnpoly(:) ltpoly(:)], [Grd.lon_rho(:) Grd.lat_rho(:)]);
    d = distance(Grd.lat_rho(:), Grd.lon_rho(:), xy(:,2), xy(:,1), ell)./1000;
    
    % Assign 1 to points in polygon, 0 to points on land, and e-folding
    % weight to all others
    
    isin = inpolygon(Grd.lon_rho, Grd.lat_rho, lnpoly, ltpoly);
 
    w = reshape(exp(-(1/Opt.efol) .* d), size(Grd.lon_rho));
    w(isin) = 1;
    w(~Grd.mask_rho) = 0;
    
    % normalize so full weight sums to 1
    
    frac(:,:,ii) = w./sum(w(:)); 
    
end

for ii = 1:nrus
    fprintf('Calculating distances: Russia %d/%d\n', ii, nrus);
    
    % Distance from each grid point to the river mouth
    
    d = distance(Grd.lat_rho(:), Grd.lon_rho(:), RusRiv.lat(ii), RusRiv.lon(ii), ell)./1000;
    
    % E-folding weight
    
    w = reshape(exp(-(1/Opt.efol) .* d), size(Grd.lon_rho));
    w(~Grd.mask_rho) = 0;
    
    % normalize so full weight sums to 1
    
    frac(:,:,ii+nak) = w./sum(w(:)); 
    
end

fprintf('Distributing across grid...\n');

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

cvt = 1000/(m2ft^3); % ft^3 -> kg, assuming fresh water

tsall = [Ts.alaska Ts.russia] * cvt; % kg/s

runoff = zeros(nxi, neta, nt);

for ii = 1:(nak+nrus)
    
    tmp = tsall(:,ii);
    tmp = bsxfun(@times, frac(:,:,ii)./Grd.area, permute(tmp, [2 3 1]));
    
    runoff = runoff + tmp;

end
runoff(isnan(runoff)) = 0; % edge rho-points (where grid cell area is undefined)

%% Export to file

Sroms = romsgeometryparams(Grd);

% File schema

fname = sprintf('runoff.kearney.efol%d.updated%s.nc', Opt.efol, datestr(today, 'yyyymm'));
Riv = bering10k_schema(Sroms, 'runoff');

% Add to history attribute

mfl = [mfilename '.m'];
if isempty(mfl)
    mfl = 'riversToRunoff.m';
end

ishis = strcmp({Riv.Attributes.Name}, 'history');
hisstr = sprintf('%s: River data added via %s\n%s', ...
    datestr(now, 'dd-mmm-yyyy HH:MM:SS'), mfl, Riv.Attributes(ishis).Value);
Riv.Attributes(ishis).Value = hisstr;

% Convert time to netcdf-style reference date and time since arrays

tdurbin = days(Ts.tbin - datetime(Sroms.tref, 'convertfrom', 'datenum')); % bin edges
tdur = (tdurbin(1:end-1) + tdurbin(2:end))./2; % midpoint

% Create file and add data

ncwriteschema(fname, Riv);
ncwrite(fname, 'runoff_time', tdur);
ncwrite(fname, 'Runoff', runoff);








