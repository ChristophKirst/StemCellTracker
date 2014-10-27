function ir = imfrmtIsRange(iSize, iFrmt, varargin)
%
% range = imfrmtRangeFromIndex(iSize, iFrmt, idx)
%
% description:
%      checks if indices form a multiplicative range
%
% input:
%      iSize   data size
%      iFormat data format
%      idx     indices
%
% output:
%      ir      true if ids form a range

if nargin == 2 
    ir = true;
elseif nargin > 2 && isnumeric(varargin{1})
    if length(varargin) == 1 % index range

        sub = imind2sub(iSize, varargin{1});
        
        n = 1;
        for i = 1:length(iFrmt)
            r = unique(sub(:,i))';
            %    if max(r) > iSize(i) || length(r) > iSize(i) || min(r) < 1
            %       error('imfrmtIsRange: range for %s out of bounds %s not in [%g,%g] !', iFrmt(i), var2char(r), 1, iSize(i))
            %    else
            n = n * length(r);
            %   end
            
        end
        
        ir = (n == size(sub,1));
        
    else
        ir = length(varargin{i}) == length(iFrmt) && all(cellfun(@isnumeric, varargin));
    end
else
    ir = true;
end