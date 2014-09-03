function q = overlapQuality(ap, varargin)
%
% q = overlapQuality(ap, varargin)
%
% description:
%    for overlaps with no signal accurate alignment is not possible
%    this routine determines the quality of the overlap statistics in the AlignmentPair ap 
%
% input:
%    ap     AlignmentPair class 
%    param  parameter struct with entries
%           .method               'grid' = use potential overlap region or 'aligned' use overlap regoind determined by ap.shift
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

thmax = getParameter(param, 'threshold.max', Inf);
thmin = getParameter(param, 'threshold.min', -Inf);
thvar = getParameter(param, 'threshold.var', 0);

wtmax = getParameter(param, 'weight.threshold.max', 1);
wtmin = getParameter(param, 'weight.threshold.min', 1);
wtvar = getParameter(param, 'weight.threshold.var', 1);

wdmax = getParameter(param, 'weight.difference.max', 0);
wdmin = getParameter(param, 'weight.difference.min', 0);
wdvar = getParameter(param, 'weight.difference.var', 0);


bmax = (thmax < Inf && wtmax ~= 0) || wdmax ~= 0;
bmin = (thmin > -Inf && wtmin ~= 0) || wdmin ~= 0;
bvar = (thvar > 0 && wtvar ~= 0) || wdvar ~= 0;

m = strcmp(getParameter(param, 'method', 'grid'), 'grid');

b = bmax || bmin || bvar;

np = length(ap);
q = zeros(1,np);

if ~b
   return
end


for p = 1:np
   
   if m
      [ov1, ov2] = ap(p).overlapOnGrid(param);
   else
      [ov1, ov2] = ap(p).overlap();
   end
   
   ov1 = double(ov1); ov2 =double(ov2);
   
   qq = 0;
   
   if bmax
      o1max = max(ov1(:)); o2max = max(ov2(:));
   
      if thmax < Inf
         if o1max < thmax
            qq = qq + wtmax * (thmax - o1max);
         end
      
         if o2max < thmax
            qq = qq + wtmax * (thmax - o2max);
         end
      end
      
      if wdmax ~= 0
         qq = qq + wdmax * abs(o1max - o2max);
      end
   end
      
   if bmin
      o1min = min(ov1(:)); o2min = min(ov2(:));
   
      if thin > -Inf
         if o1min > thmin
            qq = qq + wtmin * (o1min - thmin);
         end
      
         if o2min > thmin
            qq = qq + wtmin * (o2min - thmin);
         end
      end
      
      if wdmin ~= 0
         qq = qq + wdmin * abs(o1min - o2min);
      end
   end 
   
   if bvar
      o1var = var(ov1(:)); o2var = var(ov2(:));
   
      if thvar > 0
         if o1var < thvar
            qq = qq + wtvar * (thvar - o1var);
         end
         
         if o2var < thvar
            qq = qq + wtvar * (thvar - o2var);
         end
      end
      
      if wdvar ~= 0
         qq = qq + wdvar * abs(o1var - o2var);
      end
   end  
  
   q(p) = - qq;

end
   