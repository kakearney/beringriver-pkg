function h = plotwatersheddata(Fuw, Cst, Wbd, G, S, Nhd, R, ii, lblflag)
%PLOTWATERSHEDDATA Used to manually analyze river dataset
%
% h = plotwatersheddata(Fuw, Cst, Wbd, G, S, Nhd, R, ii, lblflag)
%
% This function is used to determine which rivers gauging stations to
% include in my Bering Sea runoff product.  It plots:
% - a zoomed out map indicating where the hydrologic unit is located
% - a set of networks showing whether certain rivers are upstream of
%   others in the HU network
% - a large map showing all HUs, the ToHUC network of those HUs, river
%   flowlines, and location of all USGS stations in watershed
% - discharge timeseries from all USGS stations in watershed
%
% All input variables come from the watersheds.m script (or the first cell
% of riversfromwatersheds.m, for Cst), except:
%
%   ii:         index of watershed to analyze
%
%   lblflag:    logical scalar, true to use smart text label positioning
%               (slower but more legible), false to use default Matlab text
%               positioning.

%-----------------------------------
% Map of where this region is
%-----------------------------------

Tmp = Wbd(Fuw.upstreamidx{ii});
latlim = minmax([Tmp.Lat]);
if diff(latlim) < 0.5
    latlim = mean(latlim) + [-0.25 0.25];
end
lonlim = minmax(wrapTo360([Tmp.Lon]));
if diff(lonlim) < 1
    lonlim = mean(lonlim) + [-0.5 0.5];
end

ltbox = latlim([1 2 2 1 1]);
lnbox = lonlim([1 1 2 2 1]);
[ltbox, lnbox] = interpm(ltbox, lnbox, 0.1);

h.ak = plotgrid('size', [1 1], 'mar', 0);
usamap('alaska');
borders('alaska', 'facecolor', rgb('light gray'), 'edgecolor', 'none');
patchm(ltbox, lnbox, 'r', 'edgecolor', 'k', 'facealpha', 0.5);
setm(h.ak.ax, 'labelformat', 'none', 'mlabellocation', 10);
  
%-----------------------------------
% HU subgraph
%-----------------------------------
   
cmap = cptcmap('Set2_08');
nts = height(Fuw.sitedata{ii});
cmaptmp = repmat(cmap, ceil(nts/8),1);

sg = subgraph(G, Fuw.upstreamidx{ii});
gidx = findnode(sg, Fuw.sitedata{ii}.huc);
hidx = findnode(sg, G.Nodes.Name(Fuw.mouthidx(ii)));

[tf,loc] = ismember(sg.Nodes.Name, {Wbd.HUC12});
sg.Nodes.Lat = arrayfun(@(X) nanmean(X.Lat), Wbd(loc)); % not perfect, but okay
sg.Nodes.Lon = arrayfun(@(X) nanmean(X.Lon), Wbd(loc));

sname = cellfun(@(a,b) sprintf('%s-%s', a, b), Fuw.sitedata{ii}.code, Fuw.sitedata{ii}.rivname, 'uni', 0);
[nidx, lbltmp] = aggregate(gidx, sname, @(x) sprintf('%s,',x{:}));
ndlbl = cell(size(sg.Nodes.Name));
[ndlbl{:}] = deal('');
ndlbl(nidx) = cellfun(@(x) x(1:end-1), lbltmp, 'uni', 0);
sg.Nodes.Label = ndlbl;

% which sites are upstream of others?

isup = false(size(gidx));

allpath = cell(0);
for ig = length(gidx):-1:1
    pth = shortestpath(sg, gidx(ig), hidx);
    
    [tf,loc] = cellfun(@(x) ismember(pth,x), allpath, 'uni', 0);
    issub = false(size(tf));
    for ip = 1:length(tf)
        issub(ip) = all(tf{ip}) & all(diff(loc{ip}) == 1);
    end
    if ~any(issub)
        allpath = [allpath; {pth}];
    end
    isup(ig) = any(ismember(pth(2:end-1), gidx));
end

nax = length(allpath);
nr = ceil(sqrt(nax));
nc = ceil(nax./nr);
h.graph = plotgrid('size', [nr nc], 'sp', 0, 'mar', 0);
for igrp = 1:length(allpath)
    g = subgraph(sg, allpath{igrp});
    
    h.graph.hg(igrp) = plot(h.graph.ax(igrp), g, 'xdata', g.Nodes.Lon, 'ydata', g.Nodes.Lat);
    labelnode(h.graph.hg(igrp), 1:numnodes(g), g.Nodes.Label);
    highlight(h.graph.hg(igrp), find(~cellfun('isempty', g.Nodes.Label))', 'nodecolor', 'r');
end
set(h.graph.ax, 'visible', 'off');%, 'xlim', lonlim, 'ylim', latlim);
expandAxes(h.graph.ax(1:nax));

% pathlbl = cellfun(@(x) [sprintf('%s->', ndlbl{x(1:end-1)}) ndlbl{x(end)}], allpath, 'uni', 0);


% h2.fig = figure('color', 'w');
% setpos(h2.fig, '# # 800px 600px');
% h2.ax = axes('position', [0 0  1 1]);
% h2.g = plot(sg, 'layout', 'force', 'edgecolor', rgb('light gray'), 'nodecolor', rgb('gray'));
% highlight(h2.g, gidx, 'nodecolor', 'r', 'markersize', 3);
% highlight(h2.g, hidx, 'nodecolor', rgb('blue'), 'markersize', 10);
% [h2.t,h2.a] = pointslabel(h2.g.XData(gidx), h2.g.YData(gidx), Fuw.sitedata{ii}.code, 'fontsize', 8);
% set(h2.t, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));
% set(h2.a, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));
% set(h2.ax, 'visible', 'off');
% 
% % annotation('textbox', [0 0 1 1], 'string', pathlbl, ...
% %     'vert', 'top', 'horiz', 'left', 'fontsize', 8);
% 
% fprintf('%s\n', pathlbl{:});

%-----------------------------------
% Map of hydrologic units and rivers
%-----------------------------------

spec = makesymbolspec('Polygon', ...
    {'Default', 'edgecolor', 'w'}, ...
    {'HUTYPE', 'S', 'facecolor', rgb('light pink')}, ...   % standard (drains to point)
    {'HUTYPE', 'C', 'facecolor', rgb('light gray')}, ...   % closed basin
    {'HUTYPE', 'F', 'facecolor', rgb('light cyan')}, ...   % frontal
    {'HUTYPE', 'M', 'facecolor', rgb('light purple')}, ... % multiple outlet
    {'HUTYPE', 'W', 'facecolor', rgb('light blue')}, ...   % water
    {'HUTYPE', 'I', 'facecolor', rgb('light tan')});       % island


[~,sloc] = ismember(Fuw.sites{ii}, S.code);

% h = plotgrid('size', [1 1], 'sp', 0, 'mar', 0.02, 'mt', 0.2);
h.map = plotgrid('size', [1 1], 'sp', 0, 'mar', 0.01, 'mt', 0);
set(h.map.fig, 'position', get(0, 'screensize')+[5 5 -10 -10]);

usamap(latlim, lonlim);
setm(h.map.ax, 'flinewidth', 1);
geoshow(Cst, 'edgecolor', 'none', 'facecolor', rgb('light gray'));
geoshow(Tmp, 'symbolspec', spec, 'facealpha', 0.6);
geoshow(Tmp(1), 'edgecolor', 'k', 'facecolor', 'none', 'linewidth', 1);
geoshow(Nhd, 'color', rgb('medium blue'));
geoshow(R, 'color', rgb('medium blue'));
plotm(S.lat(sloc), S.lon(sloc), 'marker', 'o', 'markerfacecolor', 'r', ...
    'markeredgecolor', 'none', 'markersize', 3, 'linestyle', 'none');

[xsite, ysite] = mfwdtran(S.lat(sloc), S.lon(sloc));
sname = regexp(S.name(sloc), '\s(AT|NR|BL|AB|ON|IN|[\d\.]*MI|[\d\.]* MI)\s', 'split');
sname = cellfun(@(x) x{1}, sname, 'uni', 0);
sname = cellfun(@(a,b) sprintf('%s-%s',a,b), sname, S.code(sloc), 'uni', 0);

if nargin < 9
    lblflag = true;
end
if lblflag
    [ht, ha] = pointslabel(xsite, ysite, sname, 'fontsize', 8, 'radius', 0, 'fontweight', 'b');
    set(ha, 'color', 'r');
else
   ht = text(xsite, ysite, sname,  'fontsize', 8, 'fontweight', 'b');
end


%     textm(S.lat(sloc), S.lon(sloc), S.code(sloc), 'vert', 'base', 'horiz', 'center', 'fontsize', 4);
title(Tmp(1).NAME);

setm(h.map.ax, 'fontsize', 6)

[xnd,ynd] = mfwdtran(sg.Nodes.Lat, sg.Nodes.Lon);
plot(sg, 'xdata', xnd, 'ydata', ynd, 'edgecolor', 'g', 'nodecolor', rgb('green'));

%-----------------------------------
% Timeseries
%-----------------------------------
    

isd = ~cellfun('isempty', Fuw.sitedata{ii}.discharge);
ttmp = cell(size(Fuw.sitedata{ii}.discharge));
dtmp = cell(size(ttmp));
[ttmp{:}] = deal(datetime(1950,1,1));
[dtmp{:}] = deal(NaN);

ttmp(isd) = Fuw.sitedata{ii}.timed(isd);
dtmp(isd) = Fuw.sitedata{ii}.discharge(isd);

h.ts = plotgrid('function', {@(x,y) plot(x,y/1000, '.'), ttmp, dtmp}, ...
    'sv', -0.05, 'outputs', {'ln'}, 'mar', 0.02);
set(h.ts.fig, 'color', 'w');
set(h.ts.ax(1:end-1), 'xcolor', 'none');
h.ts.ln = cat(1, h.ts.ln{:});
set(h.ts.ln, 'markersize', 1);

xlim = [datetime(1950,1,1) datetime('today')];
ylim = [0 prctile(cell2mat(dtmp)/1000, 98)];

set(h.ts.ax, 'clipping', 'off', 'box', 'off', 'color', 'none', 'box', 'off', 'tickdir', 'out', ...
    'ylim', ylim, ...
    'xlim', xlim);
set(h.ts.ax(2:2:end), 'yaxisloc', 'right');

sname = cellfun(@(a,b) sprintf('%s-%s', a, b), Fuw.sitedata{ii}.code, Fuw.sitedata{ii}.rivname, 'uni', 0);
h.ts.txt = multitextloc(h.ts.ax, sname, 'southwest');

set(h.ts.ln, {'color'}, num2cell(cmaptmp(1:nts,:),2));
set(h.ts.ax, {'ycolor'}, num2cell(cmaptmp(1:nts,:),2));
set(h.ts.txt, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));




