function [driv, filltype] = slmfitrivers(dsites, driv, rivname, sitecode1, sitecode2, mode, rivlist, dcode, filltype, varargin)
%
%   sitecode1:  downstream river site code (y coordinate in fit)
%
%   sitecode2:  upstream river site code (x coordinate in fit)

[~,loc] = ismember({sitecode1, sitecode2}, dcode);
[~,idx] = ismember(rivname, rivlist);

y = dsites(:, loc(1));
x = dsites(:, loc(2));

hasboth = ~isnan(x) & ~isnan(y);

switch mode
    case 'fit'
        slm = slmfit(x(hasboth), y(hasboth), varargin{:});
    case 'engine'
        slm = slmengine(x(hasboth), y(hasboth), varargin{:});
        
        fill = isnan(driv(:,idx)) & ~isnan(x);
        driv(fill,idx) = slmeval(x(fill), slm);
        filltype(fill,idx,1) = 2;
        filltype(fill,idx,2) = loc(2);
end

