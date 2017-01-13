function [driv, filltype] = fillriver(dsites, driv, rivname, sitecode, rivlist, dcode, filltype)
%FILLRIVER Plug data from site timeseries into gaps of river timeseries

[tf,idx] = ismember(rivname, rivlist);
[tf,loc] = ismember(sitecode, dcode);

ismissing = isnan(driv(:,idx)) & ~isnan(dsites(:,loc));

driv(ismissing,idx) = dsites(ismissing,loc);
filltype(ismissing,idx,1) = 1;
filltype(ismissing,idx,2) = loc;
