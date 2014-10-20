function rgs = imfrmtRangeFromFormatAndVarargin(frmt, varargin)
%
% rgs = imfrmtRangeFromFormatAndVarargin(frmt, param)
%
% description: 
%    parses parameter in varargin and
%    keeps those in format frmt
%    uses last format specification if multiple
%
% input:
%    frmt    format
%    param   parameter

rgs = parseParameter(varargin);

% keep indices of the format entries only
fnames = fieldnames(rgs);
ids = ismember(lower(fnames), num2cell(lower(frmt)));
rgs = rmfield(rgs, fnames(~ids));

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
end

end