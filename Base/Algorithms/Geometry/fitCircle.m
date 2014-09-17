function [c,r, varargout] = fitCircle(points)
%
% [c,r, err] = fitCircle(points)
%
% description:
%    uses algebraic expression to describe a cricle that is linear in the circle parameters 
%    and fits best least squares parameters using singular value decomposition
%
% input
%    points     array of points first row is x second is y
%
% output
%    c          center as column
%    r          radius 
%    err        error given by singular value 


% circle given by: 0 == p x.x + q.x + s
% gives linear system Bu == 0. we can constrain u = (p,q,s) to ||u|| = conts > 0 
% to solve this we find the singular right eigenvecotr with smallest eigenvalue


if size(points,1) ~= 2
   error('firCircle: input array wrong dimensions, expects coordinates along column direction');
end

x = points(1,:)';
y = points(2,:)';
n = length(x);

B = [x.^2 + y.^2, x, y, ones(n,1)];
[~,S,V] = svd(B);

se1 = V(:,end);
p = se1(1);
q = se1(2:3);
s = se1(4);
  
% center and radfrom optimal parameter
c = -q/(2*p);
r = sqrt(sum(q.^2)/(4*p^2)-s/p);

if nargout > 0
   k = size(V,2);
   varargout{1} = S(k,k);
end

end