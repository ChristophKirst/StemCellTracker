function rgs = imfrmtRangeFromVarargin(varargin)
%
% rgs = imfrmtRangeFromVarargin(param)
%
% description: 
%    parses parameter and uses last format specification if multiple
%
% input:
%    param parameter
%
% output:
%    rgs   ranges as in input. for equal entries up to upper/lower only last is kept

rgs = parseParameter(varargin);

% remove all but last coordinate specs
fnames = fieldnames(rgs);
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
   rgs = rmfield(rgs, fnames(unique(rmids)));
   %fnames = fieldnames(rgs);
end

end