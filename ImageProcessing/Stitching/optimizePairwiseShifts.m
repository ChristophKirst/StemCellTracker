function [pwdist, icenters] = optimizePairwiseShifts(pwdist)
%
% shifts = optimizePairwiseShifts(pwdist)
%
% description:
%    use least squares optimization to find globally optimal shifts from pairwise distances
%    the error function is sum (x_i + s_ij - x_j)^2
%
% input:
%    pwdist    struct array with pairwise shift imformation (.from .to .shift)
%
% output:
%    pwdistopt otimized consistant shifts as struct (.from .to .shift)
%
% See also: alignImagesOnGrid

npairs = length(pwdist);

if npairs <= 0
   icenters = [];
   return
end

% construct the mappings between node ids and index 1:nimages

indexToNodes = unique([[pwdist.from], [pwdist.to]]);

nimages = length(indexToNodes);

nodeToIndex = zeros(1, max(indexToNodes));
nodeToIndex(indexToNodes) = 1:nimages;


dim = length(pwdist(1).shift);

n = dim * npairs;
m = dim * (nimages - 1); % first image is fixed to be at [0,0]

% derivative of the error gives constraints y - X b == 0
% y are the shifts, b the centers of the images, X is derived from the error terms

%y
y = zeros(n, 1);
k = 1;
for i = 1:npairs
   sh = pwdist(i).shift;
   for d = 1:dim
      y(k) = sh(d);
      k = k + 1;
   end
end


% X
X = zeros(n,m);
k = 1;
for i = 1:npairs
   for d = 1:dim  
      pwdfromIdx = nodeToIndex(pwdist(i).from);
      if pwdfromIdx > 1
         X(k, (pwdfromIdx - 2) * dim + d) = -1;
      end
      pwdtoIdx = nodeToIndex(pwdist(i).to);
      if pwdtoIdx > 1
         X(k, (pwdtoIdx - 2) * dim + d) = 1;
      end
      k = k + 1;
   end
end

% find the centers of the images via pseudo inverse

icenters = pinv(X) * y;
icenters = [zeros(1,dim), icenters'];
icenters = reshape(icenters, dim, []);

%icenters
%size(icenters)


% translate back to consistent optimized shifts

for i = 1:npairs
   fromid = nodeToIndex(pwdist(i).from);
   toid = nodeToIndex(pwdist(i).to);
   pwdist(i).shift = (icenters(:,toid) - icenters(:,fromid));
end

end
