function sub = imind2sub(siz, idx)
%
% sub = imind2sub(siz, idx)
%
% description: 
%     returns array of position coordinates for the linear indices idx
%     assuming image of size siz
%
% input:
%     siz      size vector
%     idx      linear indices
%
% output:
%     sub      array of position indices [idim1, idim2, ...]
%
% note:
%    ind2sub only returns d vectors, here we return single matrix
%    pixel coordinates are the transpose of sub
%
% See also: ind2sub, imsub2ind

siz = double(siz);
lensiz = length(siz);
sub = zeros(length(idx), lensiz);

k = [1 cumprod(siz(1:end-1))];
for i = lensiz:-1:1,
     vi = rem(idx-1, k(i)) + 1;
     vj = (idx - vi)/k(i) + 1;
     sub(:, i) = vj;
     idx = vi;
end

end