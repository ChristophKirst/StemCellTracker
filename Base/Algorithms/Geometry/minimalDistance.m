function d = minimalDistance(ref, points)
%
% d = minimalDistance(ref, points)
%
% description:
%    finds minimal distance of rererence point to the points
%
% input
%    ref        d dim vector or d x m array of reference points
%    points     array of d dim points (d x n array)
%
% output
%    d          distances m dim vector        

%nrefs = size(ref, 2);
%npoints = size(points, 2);

d = distanceMatrix(ref, points);
d = min(d,[],2);

end