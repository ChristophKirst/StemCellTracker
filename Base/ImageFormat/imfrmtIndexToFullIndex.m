function idx = imfrmtIndexToFullIndex(iSize, iFrmt, iRange, varargin)
%
% range = imfrmtIndexFromIndex(iSize, iFrmt, iRange, varargin)
%
% description:
%      converts index w.r.t. restricted range to full index assuming size and format of iSize and iFrmt
%
% input:
%      iSize   full data size
%      iFrmt   full data format
%      iRange  range spcifying the restritions
%
% output:
%      idx     linear indices

rSize = imfrmtRangeSize(iSize, iFrmt, iRange);

% parse indices
if nargin == 3
   error('imfrmtIndexToFullIndex: no index specified!');
elseif nargin == 4
   idx = varargin{1};
   sub = imind2sub(rSize, idx);
else
   if length(varargin) ~= length(iFrmt)
      error('imfrmtIndexToFullIndex: dimension mistmatch!');
   end
   sub = ndgridc(varargin{:});
   sub = cellfunc(@(x)x(:), sub);
   sub = [sub{:}];
end

% convert sub indices via ranges to sub full indices
if ~isemptystruct(iRange)
   n = length(iFrmt);
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
            sub(:, ii) = rSize(ii) - sub(:, ii) + 1;
         end
            
         if lower(rf) == rf
            r = iSize(ii) - r + 1;
         end
         
         % convert indices
         sub(:, ii) = r(sub(:, ii));
         
         % back to output format
         if lower(f) == f
            sub(:, ii) = iSize(ii) - sub(:, ii) + 1;
         end
      end
   end
end

% convert to full linear indices
idx = imsub2ind(iSize, sub)';

end