function range = imfrmtRangeFromIndex(iSize, iFrmt, idx)
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
%
% See also: tagExpression

range = struct();
sub = imind2sub(iSize, idx);

n = 1;
for i = 1:length(iFrmt)
   r = unique(sub(:,i))';
   if max(r) > iSize(i) || length(r) > iSize(i)
      error('imfrmtRangeFromIndex: inconsistent ranges !')
   elseif length(r) == iSize(i)
      n = n * length(r);
   else   
      range.(iFrmt(i)) = r;
      n = n * length(r);
   end

end   

if n ~= length(idx)
   error('imfrmtRangeFromIndex: ranges not multiplicative !')
end