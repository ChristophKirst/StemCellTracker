function pol = polygonDilate(pol, dist)
%
% pol = polygonDilate(pol, dist)
%
% description:
%    increase a polygon by adding a border of size bwidth
%
% input:
%    pol     array of polygon coordinates as row vectors
%    dist    width of the border to dilate / erode
%
% ouput:
%    pol     polygon with added border
%
% See also: bufferPolygon

pol = bufferPolygon(pol, dist);

if length(pol) == 1
   pol = pol{1};
end

end