function a = triangulationArea(tri)
%
% a = triangulationArea(pol)
%
% description:
%   calculatate area of a triangulation (this should be mathworks job really!)
%
% input:
%   tri    triangulation object
%
% output:
%   a      area
%
% See also: polygonArea

T = tri.ConnectivityList;
P = tri.Points;

T1 = T(:,1); T2 = T(:,2); T3 = T(:,3);
X = P(:,1); Y = P(:,2);

x1 = X(T1); x2 = X(T2); x3 = X(T3);
y1 = Y(T1); y2 = Y(T2); y3 = Y(T3);

a = sum(1/2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)));

end