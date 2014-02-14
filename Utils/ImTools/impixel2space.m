function xyz = impixel2space(isize, extend, hwl)
%
% xyz = impixel2space(iszie, extend, hwl)
%
% description:
%     converts pixel coordinates to spatial coordinates (no xy exchange) 
%
% input:
%     isize     image pixel size [h, w, l]
%     extend    spatial extends [xmin, ymin, zmin; xmax, ymax, zmax]
%     hwl       pixel coordinates to convert in the form [h; w; l]
%
% output:
%     xyz       spatial coordinates [x; y; z]
%
% See also: imspace2pixel


ssize = extend(2,:) - extend(1,:);
fac = ssize ./ isize;
xyzmin = extend(1,:);
np = size(hwl,1);

xyz = (hwl - 0.5) .* repmat(fac, np, 1) + repmat(xyzmin, np,1);

end