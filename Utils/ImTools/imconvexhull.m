function [xyz,tri] = imconvexhull(labelimage)
%
% [xyz,tri] = imconvexhull(labelimage)
%
% description:
%     calculates convexhull for each label in labeld image labelimage in pixel coordinates
% 
% input:
%     labelimage   the labeled image
%
% output:
%     xyz,tri      vertices and faces of the convex hulls as cell aarrays
%
% See also: bwconvexhull3d

labels = imlabel(labelimage);
n = length(labels);

xyz = cell(n);
tri = cell(n);

i = 1;
for l = labels
    disp(l)
    bw = labelimage == l;
    [t, x] = bwconvhull3d(bw);
    tri{i} = t; xyz{i} = x;
    i = i + 1;
end
    
end
    

