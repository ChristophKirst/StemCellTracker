function range = imfrmtRangeFromRangeAndVarargin(refRange, varargin)
%
% range = imfrmtRangeFromRangeAndVarargin(range, param)
%
% description: 
%    overwrites ranges in range by the ranges specifed in varargin and drops all others varargins
%    uses last format specification if multiple
%
% input:
%    refRange   reference range
%    param      parameter
%
% output:
%    range      parsed range

range = parseParameter(varargin);

% keep indices of the format entries only

refNames = fieldnames(refRange);
fNames   = fieldnames(range);

ids = ismember(lower(fNames), lower(refNames));
range = rmfield(range, fNames(~ids));

% remove all but last coordinate specs
fnames = fieldnames(range);
fnamesl = lower(fnames);
rmids = [];
for i = 1:length(fnames)
   id = ismember(fnamesl, fnamesl{i});
   id = find(id);
   if length(id) > 1
      rmids = [rmids, id(1:end-1)]; %#ok<AGROW>
   end
end
if ~isempty(rmids)
   range = rmfield(range, fnames(unique(rmids)));
   %fnames = fieldnames(range);
end

end