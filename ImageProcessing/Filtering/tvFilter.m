function img = tvFilter(img, varargin)
%
% img = tvFilter(img, param)
%
% description:
%   Split Bregman Anisotropic Total Variation Denoising
%   minimizing  1/2||img-imgf||_2^2 + mu*ATV(imgf)
% 
% input:
%   img   image
%   param parameter struct with entries
%         mu     weight of TV w.r.t. mean error
%
% Refs:
%  *Goldstein and Osher, The split Bregman method for L1 regularized problems
%   SIAM Journal on Imaging Sciences 2(2) 2009
%  *Micchelli et al, Proximity algorithms for image models: denoising
%   Inverse Problems 27(4) 2011

g = img(:);
n = length(g);
[B Bt BtB] = DiffOper(sqrt(n));
b = zeros(2*n,1);
d = b;
u = g;
err = 1;k = 1;
tol = 1e-3;
lambda = 1;
while err > tol
    fprintf('it. %g ',k);
    up = u;
    [u,~] = cgs(speye(n)+BtB, g-lambda*Bt*(b-d),1e-5,100); 
    Bub = B*u+b;
    d = max(abs(Bub)-mu/lambda,0).*sign(Bub);
    b = Bub-d;
    err = norm(up-u)/norm(u);
    %fprintf('err=%g \n',err);
    k = k+1;
end
%fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);

img = reshape(u, size(img));

end

function [B Bt BtB] = DiffOper(N)
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;
end