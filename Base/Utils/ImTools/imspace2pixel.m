function pql = imspace2pixel(isize, extend, xyz)
%
% xyz = impixel2space(iszie, extend, pql)
%
% description:
%     converts pixel coordinates to spatial coordinates (no xy exchange) 
%
% input:
%     isize     image pixel size [p, q, l]
%     extend    spatial extends [xmin, ymin, zmin; xmax, ymax, zmax]
%     xyz       spatial coordinates [x; y; z]
%
% output:
%     pql       pixel coordinates to convert in the form [p; q; l]
%
% See also: imspixel2space

ssize = extend(1,:) - extend(2,:);
fac = isize ./ ssize;
xyzmin = extend(2,:);
np = size(xyz,1);

pql = round((xyz -repmat(xyzmin, np,1)).* repmat(fac, np, 1)) + 0.5;

end
