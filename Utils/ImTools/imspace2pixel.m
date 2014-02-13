function hwl = imspace2pixel(isize, extend, xyz)
%
% xyz = impixel2space(iszie, extend, hwl)
%
% description:
%     converts pixel coordinates to spatial coordinates (no xy exchange) 
%
% input:
%     isize     image pixel size [h, w, l]
%     extend    spatial extends [xmin, ymin, zmin; xmax, ymax, zmax]
%     xyz       spatial coordinates [x; y; z]
%
% output:
%     hwl       pixel coordinates to convert in the form [h; w; l]
%
% See also: imspixel2space

ssize = extend(1,:) - extend(2,:);
fac = isize ./ ssize;
xyzmin = extend(2,:);
np = size(xyz,1);

hwl = round((xyz -repmat(xyzmin, np,1)).* repmat(fac, np, 1));

end
