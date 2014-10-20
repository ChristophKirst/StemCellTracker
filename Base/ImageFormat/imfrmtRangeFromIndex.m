function range = imfrmtRangeFromIndex(iSize, iFrmt, varargin)
%
% range = imfrmtRangeFromIndex(iSize, iFrmt, idx)
%
% description:
%      converts index to a range
%
% input:
%      iSize   data size
%      iFormat data format
%      idx     indices
%      param   (optional) parameter struct with entries
%              .check    check for multiplicative tag grid (false)
%
% output:
%      range   struct with .name = range entries


range = struct();
if nargin == 3
   sub = imind2sub(iSize, varargin{1});
else
   sub = ndgridc(varargin{:});
   sub = cellfunc(@(x)x(:), sub);
   sub = [sub{:}];
end

n = 1;
for i = 1:length(iFrmt)
   r = unique(sub(:,i))';
   if max(r) > iSize(i) || length(r) > iSize(i)
      error('imfrmtRangeFromIndex: range for %s out of bounds %s not in [%g,%g] !', iFrmt(i), var2char(r), 1, iSize(i))
   elseif length(r) == iSize(i)
      n = n * length(r);
   else   
      range.(iFrmt(i)) = r;
      n = n * length(r);
   end

end   

if n ~= size(sub,1)
   error('imfrmtRangeFromIndex: ranges not multiplicative !')
end