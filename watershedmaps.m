% This script plots the watersheds and all USGS sites located in each 

%% Coastlines

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
for ir = 1:5
    rfile = fullfile(gshhsdir, 'WDBII_shp', lev, sprintf('WDBII_river_%s_L%02d.shp', lev, ir));
    Riv{ir} = shapeprjread(rfile, 'BoundingBox', [aklonlim' aklatlim']);
end

%% Plots of each watershed 

spec = makesymbolspec('Polygon', ...
    {'Default', 'edgecolor', 'w'}, ...
    {'HUTYPE', 'S', 'facecolor', rgb('light pink')}, ...   % standard (drains to point)
    {'HUTYPE', 'C', 'facecolor', rgb('light gray')}, ...   % closed basin
    {'HUTYPE', 'F', 'facecolor', rgb('light cyan')}, ...   % frontal
    {'HUTYPE', 'M', 'facecolor', rgb('light purple')}, ... % multiple outlet
    {'HUTYPE', 'W', 'facecolor', rgb('light blue')}, ...   % water
    {'HUTYPE', 'I', 'facecolor', rgb('light tan')});       % island

if ~exist('watershedplots', 'dir')
    mkdir('watershedplots');
end
for ii = 1:height(Fuw)
    Tmp = Wbd(Fuw.upstreamidx{ii});
    latlim = minmax([Tmp.Lat], 'expand');
    if diff(latlim) < 0.5
        latlim = mean(latlim) + [-0.25 0.25];
    end
    lonlim = minmax([Tmp.Lon], 'expand');
    if diff(lonlim) < 1
        lonlim = mean(lonlim) + [-0.5 0.5];
    end
    
    [~,sloc] = ismember(Fuw.sites{ii}, S.code);
    
    h = plotgrid('size', [1 1], 'sp', 0, 'mar', 0.02, 'mt', 0.2);
    h.ax(2) = axes('position', [0.02 0.8 0.2 0.2]);
    h.ax = h.ax([2 1]);
    
    setpos(h.fig, '# # 800 600');
    set(h.fig, 'color', 'w');
    axes(h.ax(1));
    usamap('alaska');
    borders('alaska', 'facecolor', rgb('light gray'), 'edgecolor', 'none');
    patchm(latlim([1 2 2 1 1]), lonlim([1 1 2 2 1]), 'w', 'edgecolor', 'k', 'facealpha', 0.5);
    setm(h.ax(1), 'labelformat', 'none', 'mlabellocation', 10);
    
    axes(h.ax(2));
    usamap(latlim, lonlim);
    setm(h.ax(2), 'flinewidth', 1);
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
    
    [ht, ha] = pointslabel(xsite, ysite, sname, 'fontsize', 6, 'radius', 0, 'fontweight', 'b');
    set(ha, 'color', 'r');
    
%     textm(S.lat(sloc), S.lon(sloc), S.code(sloc), 'vert', 'base', 'horiz', 'center', 'fontsize', 4);
    title(Tmp(1).NAME);
    
    setm(h.ax(1), 'fontsize', 6, 'frame', 'off')
    setm(h.ax(2), 'fontsize', 6)
    
    imgfile = fullfile('watershedplots', sprintf('fuw%02d',ii));
    export_fig(imgfile, '-png', '-r300', '-nocrop')
    
    close(h.fig);
    
end
    
%% Choose-the-station plots
%
% These plots show the various timeseries, and the graph network for where
% these sites fall in the HU hierarchy.  These are used to manually decide
% which stations to analyze.

cmap = cptcmap('Set2_08');

for ii = 1%  1:height(Fuw)
    
    % Timeseries
    
    h = plotgrid('function', {@(x,y) plot(x,y/1000, '.'), Fuw.sitedata{ii}.timed, Fuw.sitedata{ii}.discharge}, ...
        'sv', -0.05, 'outputs', {'ln'});
    set(h.fig, 'color', 'w');
    set(h.ax(1:end-1), 'xcolor', 'none');
    h.ln = cat(1, h.ln{:});
    nts = length(h.ln);
    set(h.ln, 'markersize', 1);
    
    set(h.ax, 'clipping', 'off', 'box', 'off', 'color', 'none', 'box', 'off', 'tickdir', 'out', ...
        'ylim', [0 prctile(cell2mat(Fuw.sitedata{ii}.discharge)/1000, 98)], ...
        'xlim', [datetime(1950,1,1) datetime('today')]);
    set(h.ax(2:2:end), 'yaxisloc', 'right');
    
    sname = cellfun(@(a,b) sprintf('%s-%s', a, b), Fuw.sitedata{ii}.code, Fuw.sitedata{ii}.rivname, 'uni', 0);
    h.txt = multitextloc(h.ax, sname, 'southwest');
    
    cmaptmp = repmat(cmap, ceil(nts/8),1);
    set(h.ln, {'color'}, num2cell(cmaptmp(1:nts,:),2));
    set(h.ax, {'ycolor'}, num2cell(cmaptmp(1:nts,:),2));
    set(h.txt, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));
    
    % Subgraph
   
    sg = subgraph(G, Fuw.upstreamidx{ii});
    gidx = findnode(sg, Fuw.sitedata{ii}.huc);
    hidx = findnode(sg, G.Nodes.Name(Fuw.mouthidx(ii)));
    
    h2.fig = figure('color', 'w');
    setpos(h2.fig, '# # 800px 600px');
    h2.ax = axes('position', [0 0  1 1]);
    h2.g = plot(sg, 'layout', 'force', 'edgecolor', rgb('light gray'), 'nodecolor', rgb('gray'));
    highlight(h2.g, gidx, 'nodecolor', 'r', 'markersize', 3);
    highlight(h2.g, hidx, 'nodecolor', rgb('blue'), 'markersize', 3);
    [h2.t,h2.a] = pointslabel(h2.g.XData(gidx), h2.g.YData(gidx), Fuw.sitedata{ii}.code, 'fontsize', 8);
    set(h2.t, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));
    set(h2.a, {'color'}, num2cell(brighten(cmaptmp(1:nts,:), -0.6),2));
    set(h2.ax, 'visible', 'off');
    
%     pause;
%     close(h.fig);
%     close(h2.fig);
    
end

% Results
% 


