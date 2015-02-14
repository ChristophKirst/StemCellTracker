function [B,A] = polygonToBoundaries(pol, varargin)
%
% tri = polygonToTriangulation(poly)
% 
% description: 
%      converts a polygon into a triangulation
%
% input:
%      pol   cell array of countours not self-intersecting
%
% output:
%      tri   triangulation object


% get all ccordinates and boundary edges assuming the individual contours are not intersecting and not self intersecting

param = parseParameter(varargin);
a = getParameter(param, 'all', false);

if ~iscell(pol)
   pol = {pol};
end

for i = 1:length(pol)
   if size(pol{i},1) ~= 2
      error('polygonToTriangulation: polygon coordinate size %d is not 2 in contour %d!', size(pol{i},1), i);
   end
end

xy = cell2mat(pol);
%size(pol)
%size(xy)

bd = [];
id = 1;
for i = 1:length(pol)
   s = id;
   e = id + size(pol{i}, 2) - 2;
   bdi = [(s:e)' ((s+1):(e+1))'; (e+1), s];
   bd = [bd; bdi]; %#ok<AGROW>
   id = id + size(pol{i}, 2);
end

%bd
tri = delaunayTriangulation(xy', bd);

% remove holes
if ~a
   isi = tri.isInterior;
   if all(~isi)
      tri = [];
   else
      tri = triangulation(tri.ConnectivityList(isi,:), tri.Points);
   end
end

end

