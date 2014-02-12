function [X,K] = imlabel2surface(labelimage)
%
% [X,K] = imlabel2surface(labelimage)
%
% description:
%     calculates convexhull for each label in labeld image labelimage.
% 
% input:
%     labelimage   the labeled image
%
% output:
%     X,K      cells for the convexhull data
%
% See also: convexhull


labels = unique(labelimage(:))';
if labels(1) == 0
    labels = labels(2:end);
end
n = length(labels);

X = cell(n);
K = cell(n);

i = 1;
for l = labels
    disp(l)
    bw = labelimage == l;
    [k, x] = bwconvhull3d(bw);
    K{i} = k; X{i} = x;
    i = i + 1;
end
    

    

