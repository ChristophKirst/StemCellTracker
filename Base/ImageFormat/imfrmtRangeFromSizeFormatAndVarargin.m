function rgs = imfrmtRangeFromSizeFormatAndVarargin(si, frmt, varargin)
%
% rgs = imfrmtRangeFromSizeFormatAndVarargin(si, frmt, param)
%
% description: 
%    parses parameter in varargin keeps and reformats paramter as in format frmt / size si
%    uses last format specification if multiple
%
% input:
%    si      size
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
   fnames = fieldnames(rgs);
end


% flip if required
[ids, pos] = ismember(lower(fnames), num2cell(lower(frmt)));
ids = find(ids);
pos = pos(ids);

for i = 1:length(ids)
   if fnames{ids(i)} ~= frmt(pos(i))
      v = rgs.(fnames{ids(i)});
      rgs = rmfield(rgs, fnames{ids(i)});
      rgs.(frmt(pos(ids(i)))) = si(pos(i)) - v + 1;
   end
end

end