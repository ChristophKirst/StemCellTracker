function i = imfrmtInfoFromCell(d, i)
%
% i = imfrmtInfoFromCell(d)
%
% description:
%    infers image information from a numerical array and returns it as ImageInfo class

if nargin == 1
   i = ImageInfo();
end

if ~iscell(d)
   d = {d};
end

nd = ndims1(d);
csi = size(d);
csi = csi(1:nd);
i.irawcellsize = csi;
cfrmt = 'SUVW';
i.irawcellformat = cfrmt(1:nd);

i.irawdatasize   = size(d{1});
i.irawdataformat = imfrmtFormatFromSize(i.irawdatasize);
i.irawdataclass  = class(d{1});

i.initializeCellDataSizeAndFormatFromRaw();

cl = 1;
id = find(lower(i.dataFormat) == 'c', 1, 'first');
if ~isempty(id)
   cl = i.isize(id);
end
if cl == 1
   i.icolor = {'gray'};
else
   i.icolor  = imcolorlist(cl);
end



end