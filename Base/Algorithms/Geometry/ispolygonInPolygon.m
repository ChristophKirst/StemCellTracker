function b = ispolygonInPolygon(poly1, poly2)
%
% b = ispolygonInPolygon(poly1, poly2)
%
% description:
%    determines if polygon poly1 is in polygon poly2
%
% input:
%   poly1,2  polygons as arrays of points as row vectors or ROIPolygons
%
% output:
%   b        true of poly2 is in poly1

if ~isnumeric(poly1)
   poly1 = poly1.toArray;
end

if ~isnumeric(poly2)
   poly2 = poly2.toArray;
end  

xcheck = poly1(1,:);
ycheck = poly1(2,:);

xpoly  = poly2(1,:);
ypoly  = poly2(2,:);

b = all(inpolygon(xcheck,ycheck,xpoly,ypoly));

end