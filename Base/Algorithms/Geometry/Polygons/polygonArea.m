function a = polygonArea(pol)
%
% a = polygonArea(pol)
%
% description:
%   calculatate area of polygon
%
% input:
%   pol    polygon
%
% output:
%   a      area
%
% See also: triangulationArea

a = triangulationArea(polygonToTriangulation(pol));

end