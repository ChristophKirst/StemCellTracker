function idx = imfrmtIndexFromSizeAndVarargin(iSize, varargin)
%
% range = imfrmtIndexFromIndex(iSize, iFrmt, idx)
%
% description:
%      converts index to index parsing multiple sub index ranges
%
% input:
%      iSize   data size
%
% output:
%      idx     linear indices

if nargin == 1
   idx = 1:iSize;
elseif nargin == 2
   idx = varargin{1};
else
   sub = ndgridc(varargin{:});
   sub = cellfunc(@(x)x(:), sub);
   sub = [sub{:}];
   idx = imsub2ind(iSize, sub);
end

end