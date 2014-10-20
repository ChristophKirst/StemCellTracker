function i = imfrmtInfoFromData(d, i)
%
% i = imfrmtInfoFromData(d)
%
% description:
%    infers image infomration form a nuemrical array and retunrs it as ImageInfo class

if nargin == 1
   i = ImageInfo();
end

i.irawcellsize = 1;
i.irawcellformat = 'S';

i.irawdatasize   = size(d);
i.irawdataformat = imfrmtFormatFromSize(i.irawdatasize);
i.irawdataclass  = class(d);

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