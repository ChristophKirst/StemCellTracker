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

if isempty(pol)
   a = 0;
   return
end

a = triangulationArea(polygonToTriangulation(pol));

end