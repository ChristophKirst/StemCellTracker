function u = unpadarray(a, padsize, direction)
%
% u = unpadarray(a, padsize)
%
% description:
%   removes padding specified by padsize, i.e. reverses padarray
%
% input:
%   a          array to unpad
%   padsize    size specification as in padarray
%   direction  direction, either of 'both', 'post', 'pre'
%
% output:
%   u        unpadded array
%
% See also: padarray

if nargin < 3
   direction = 'both';
end

numdims = ndims(a);
padsize = padsize(:);
if (numel(padsize) < numdims)
   padsize(ndims(a)) = 0;
end

idx   = cell(1,numdims);
for k = 1:numdims
    M = size(a,k);
    switch direction
        case 'pre'
            idx{k}   = (padsize(k)+1:M); 
        case 'post'
            idx{k}   = 1:(M-padsize(k));          
        case 'both'
            idx{k}   = (padsize(k)+1:M-padsize(k));
    end
end

u = a(idx{:});

end




