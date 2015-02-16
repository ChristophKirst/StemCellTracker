function pol = polygonFromTriangulation(tri)
%
% tri = polygonFromTriangulation(tri)
% 
% description: 
%      converts a triangulation into a polygon
%
% input:
%      tri  triangulation object
%
% output:
%      pol  polygon as cell array fo contours (even odd)


T = num2cell(tri.ConnectivityList, 2);
P = tri.Points;

pol =  cellfunc(@(x) P(x, :)', T);

pol = polygonSimplify(pol);

end

