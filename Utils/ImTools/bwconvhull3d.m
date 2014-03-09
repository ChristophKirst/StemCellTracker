function [hull, xyz] = bwconvhull3d( mask )
%
% [hull, xyz] = bwconvhull3d( mask )
%
% description:
%    calculates the convex hull of mask as bw mask 
%
% input:
%    mask   3d bw image of points
%
% output:
%    hull   convex hull as mask
%    xyz    pixel coordinates of hull
%
% See also: bwconvhull convhulln


ind = find(mask);
siz = size(mask);

[x,y,z] = ind2sub(siz, ind);
xyz = [x,y,z];

k = convhulln(xyz);

% linear constrains A x + b <= 0 
c = mean(xyz(unique(k),:));
xyz = xyz - repmat(c,[size(xyz,1) 1]);
A = NaN*zeros(size(k,1),size(xyz,2));
rc=0;
for ix = 1:size(k,1)
    F = xyz(k(ix,:),:);
    if rank(F,1e-5) == size(F,1)
        rc=rc+1;
        A(rc,:)=F\ones(size(F,1),1);
    end
end

A = A(1:rc,:);
b = ones(size(A,1),1);
b = b + A*c';

[~,I]=unique(num2str([A b],6),'rows'); % eliminate dumplicate constraints:
A = A(I,:);
b = b(I);

% get pixel values
[gy, gx, gz] = meshgrid(1:siz(1), 1:siz(2), 1:siz(3));
r = [gx(:)' ; gy(:)' ; gz(:)'];

id = all(A * r - repmat(b,1, length(r))<= 0.5);

hull = mask;
hull(id) = 1;

end

