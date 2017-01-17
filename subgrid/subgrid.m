function newax = subgrid(ax, rfrac, cfrac)
%SUBGRID Divide axis in subaxes
%
% newax = subgrid(ax, nr, nc)
% newax = subgrid(ax, rfrac, cfrac)
%
% Input variables: 
%
%   ax:     array of axes handles, each will be sub-divided as indicated
%
%   nr:     number of rows to split each axis into
%
%   nc:     number of columns
%
%   rfrac:  vector of values summing to <=1, indicating fraction of height
%           used by each new axis (bottom to top)
%
%   cfrac:  vector of values summing to <=1, indicating fraction of width
%           used by each new axis (left to right)
%
% Output variables:
%
%   newax:  nr x nc cell array, each holding an array of handles the same
%           size as ax.  Arrays are arranged to match visual arrangement of
%           axes, with newax{1,1} indicating the top left corner subaxis
%           within each old axis

% Copyright 2012 Kelly Kearney


if isscalar(rfrac) && floor(rfrac)==rfrac
    height = ones(1,rfrac)./rfrac;
elseif all(rfrac > 0 & rfrac <= 1) & sum(rfrac) <= 1
    height = rfrac;
else
    error('Input 2 must be integer or vector of fractions summing to 1');
end

if isscalar(cfrac) && floor(cfrac)==cfrac
    width = ones(1,cfrac)./cfrac;
elseif all(cfrac > 0 & cfrac <= 1) & sum(cfrac) <= 1
    width = cfrac;
else
    error('Input 3 must be integer or vector of fractions summing to 1');
end

left = [0 cumsum(width(1:end-1))];
bott = [0 cumsum(height(1:end-1))];

nrow = length(height);
ncol = length(width);

newax = cell(length(height), length(width));
if verLessThan('matlab', '8.4.0')
    [newax{:}] = deal(zeros(size(ax)));
else
    [newax{:}] = deal(gobjects(size(ax)));
end

for iax = 1:numel(ax)
  
    axpos = get(ax(iax), 'position');
    hfig = get(ax(iax), 'parent');
    
    l = axpos(1) + left.*axpos(3);
    b = axpos(2) + bott.*axpos(4);    
    b = b(end:-1:1);
    
    w = axpos(3).*width;
    h = axpos(4).*height;
    h = h(end:-1:1);
    
    for ir = 1:nrow
        for ic = 1:ncol
            newax{ir,ic}(iax) = axes('position', [l(ic) b(ir) w(ic) h(ir)], 'parent', hfig);
        end
    end
    
    
%     w = axpos(3)/sz(2);
%     h = axpos(4)/sz(1);

%     for inew = 1:numel(left)
%         newax{inew}(iax) = axes('position', [l(inew) b(inew) w(inew) h(inew)]);
%     end
    delete(ax(iax));
    
end




    
   