function isize = imfrmtRangeSize(isize, ifrmt, varargin)
%
% trg= imfrmtRangeSize(isize, ifrmt, trg)
%
% description: 
%     detemines data size given the full image size isize and the ranges rgs
%
% input:
%     isize      image size 
%     ifrmt      the reference image format used to detect entries to delete w.r.t to outfrmt
%     rgs        ranges
%
% output:
%     isize      image size given the rgs
% 
% See also: imfrmtReformatRange

rgs = parseParameter(varargin);

fnames = fieldnames(rgs);
fnamesl = lower(fnames);
ifrmtl = lower(ifrmt);

[id, pos] = ismember(fnamesl, num2cell(ifrmtl));

for i = 1:length(id)
   if id(i)
      isize(pos(i)) = length(rgs.(fnames{i}));
   end
end

end
   