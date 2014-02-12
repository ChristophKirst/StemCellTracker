function dist = distanceMatrix(r0, r1)
%
%   dist = distanceMatrix(r0, r1)
%
% description:
%   calculates pairwise distances between coordinates in arrays r0 and r1
%
% input:
%   r*        arrays of coordinates
%
% output:
%   dist      the pairwise Euclidean distances as a matrix
%
% todo:
%             compiled / parallel version 
%             for r0 == r1 speed up due to symmetry
%

if nargin < 2
   r1 = r0;
end

if isentry(r0, 'r')
   r0 = r0.r;
end
if isentry(r1, 'r')
   r1 = r1.r;
end

n0 = size(r0,2);
n1 = size(r1,2);
dim = size(r0,1);

dist = sqrt(sum(bsxfun(@minus,reshape(r0',[n0,1,dim]),reshape(r1',[1,n1,dim])).^2,3));

end
