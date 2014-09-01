function q = overlapQuality(stats, varargin)
%
% q = overlapQuality(stats, varargin)
%
% description:
%    for overlaps with no signal accurate alignment is not possible
%    this routine determines the quality of the overlap statistics stats 
%
% input:
%    stats  overlap statistics as obtained form overlapStatisticsImagePair
%    param  parameter struct with entries
%           .threshold.max        if max is below this value decrease quality  (Inf)
%           .threshold.min        if min is above this value decrease quality  (-Inf)
%           .threshold.var        if variance is below this value decrease quality (0)
%           .weight.threshold.max weights (1)
%           .weight.threshold.min (1)
%           .weight.threshold.var (1)
%           .weight.difference.max
%           .weight.difference.min
%           .weight.difference.var
%
% output:
%    quality   scalar determining the quality of the overlap
%

param = parseParameter(varargin{:});

q = 0;
ch = stats.from;
q = q + checkThresholds(ch, param);

ch = stats.to;
q = q + checkThresholds(ch, param);

w = getParameter(param, 'weight.difference.max', 0);
if w > 0
   q  = q + w * abs(stats.from.max - stats.to.max);
end

w = getParameter(param, 'weight.difference.min', 0);
if w > 0
   q  = q + w * abs(stats.from.min - stats.to.min);
end

w = getParameter(param, 'weight.difference.var', 0);
if w > 0
   q  = q + w * abs(stats.from.mvar - stats.to.var);
end

q = -q;

end


% helper
function q = checkThresholds(ch, param)
   q = 0;
   th = getParameter(param, 'threshold.max', Inf);
   if th < Inf && ch.max < th
      q = q + getParameter(param, 'weight.threshold.max', 1) * (th - ch.max);
   end
   
   th = getParameter(param, 'threshold.min', -Inf);
   if th > -Inf && ch.min > th
      q = q + getParameter(param, 'weight.threshold.min', 1) * (ch.min - th);
   end
   
   th = getParameter(param, 'threshold.var', Inf);
   if th < Inf && ch.var < th
      q = q + getParameter(param, 'weight.threshold.var', 1) * (ch.var - th);
   end
end
   