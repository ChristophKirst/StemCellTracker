function [R, T, C] = optimalTransformation(X, Y)
%
% [R, T, C] = optimalTransformation(X, Y)
%
% description: 
%    Finds the rotation matrix R, translation T and scaling C 
%    that minimize (C * R * X + T - Y).^2.
%
% input:
%    X,Y  coordinate arrays, each point as column
%
% output:
%    R    optimal rotation matrix
%    T    optimal translation
%    C    optimal scaling
%
% reference: 
%    Least-Squares Est. of Transformation Parameter Between Two Point Patterns
%    S. Umeyama, IEEE Trans Pattern Analysis Machine Inteligence 13, p376 1991
%

n = size(X,2);
dim = size(X,1);

mX = mean(X,2);
mY = mean(Y,2);

K = diag(ones(n,1)) - 1/n * ones(n,n);

vX = 1/n * sum(sum((X*K).^2));
Sigma = 1/n * Y * K * X';

[U,D,V] = svd(Sigma);

S = diag(ones(dim,1));
if det(Sigma) < 0
   S(dim,dim) = -1;
end

R = U * S * V';
C = 1/vX * trace(D * S);
T = mY - (C * R * mX);

end



