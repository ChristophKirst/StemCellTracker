function idx = imfrmtFullIndexToIndex(iSize, iFrmt, iRange, varargin)
%
% range = imfrmtFullIndexToIndex(iSize, iFrmt, iRange, varargin)
%
% description:
%      converts index w.r.t. full range assuming size and format of iSize and iFrmt 
%      to index restricted range using iRange
%
% input:
%      iSize   full data size
%      iFrmt   full data format
%      iRange  range specifying the restritions
%
% output:
%      idx     restricted linear indices


% parse indices
if nargin == 3
   error('imfrmtFullIndexToIndex: no index specified!');
elseif nargin == 4
   idx = varargin{1};
   sub = imind2sub(iSize, idx(:));
else
   if length(varargin) ~= length(iFrmt)
      error('imfrmtFullIndexToIndex: dimension mistmatch!');
   end
   sub = ndgridc(varargin{:});
   sub = cellfunc(@(x)x(:), sub);
   sub = [sub{:}];
end

% convert full sub indices via ranges to resticted sub indices
if ~isemptystruct(iRange)
   
   rSize = imfrmtRangeSize(iSize, iFrmt, iRange);
   rnames = fieldnames(iRange);
   [ids, pos] = ismember(num2cell(lower(iFrmt)), lower(rnames));
   
   if any(ids)
      ids = find(ids);
      for i = 1:length(ids)
         ii = ids(i);
         
         rf = rnames{pos(ii)};
         f = iFrmt(ii);
         r = iRange.(rf);
         
         % convert to positiv axis
         if lower(f) == f
            sub(:, ii) = iSize(ii) - sub(:, ii) + 1;
         end
            
         if lower(rf) == rf
            r = iSize(ii) - r + 1;
         end
         
         % convert indices
         
         [~, rpos] = ismember(sub(:,ii), r);
         
         if any(rpos == 0)
            error('imfrmtFullIndexToIndex: index in format %s out of range!', f);
         end
         
         sub(:, ii) = rpos;
         
         % back to output format
         if lower(f) == f
            sub(:, ii) = rSize(ii) - sub(:, ii) + 1;
         end
      end
   end
   
   % convert to full linear indices
   
   si = size(idx);
   idx = imsub2ind(rSize, sub)';   
   idx = reshape(idx, si);
   
else
 
   si = size(idx);
   idx = imsub2ind(iSize, sub)'; 
   idx = reshape(idx, si);
   
end