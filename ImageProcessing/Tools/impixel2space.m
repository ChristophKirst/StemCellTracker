function xyz = impixel2space(isize, extend, pql, scale)
%
% xyz = impixel2space(iszie, extend, pql)
%
% description:
%     converts pixel coordinates to spatial coordinates (no xy exchange) 
%
% input:
%     isize     image pixel size [p, q, l]
%     extend    spatial extends [xmin, ymin, zmin; xmax, ymax, zmax]
%     pql       pixel coordinates to convert in the form [p; q; l]
%     scale     (optional) optional scaling
%
% output:
%     xyz       spatial coordinates [x; y; z]
%
% See also: imspace2pixel

if nargin < 4
   scale =ones(1,3);
end

ssize = extend(2,:) - extend(1,:);
fac = ssize ./ isize .* scale;
xyzmin = extend(1,:);
np = size(pql,1);

xyz = (pql - 0.5) .* repmat(fac, np, 1) + repmat(xyzmin, np,1);

end