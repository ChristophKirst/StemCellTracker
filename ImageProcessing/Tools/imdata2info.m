function i = imdata2info(d)
%
% i = imdata2info(d)
%
% description:
%    infers image infomration form a nuemrical array and retunrs it as ImageInfo class

i = ImageInfo();
i.idatasize   = size(d);
i.idataformat = imsize2format(i.idatasize);
i.idataclass  = class(d);
i.irawformat  = i.idataformat;
i.irawsize    = i.idatasize;

cl = 1;
id = find(i.idataformat == 'c', 1, 'first');
if ~isempty(id)
   cl = i.isize(id);
end
if cl == 1
   i.icolor = {'gray'};
else
   i.icolor  = imcolorlist(cl);
end
i.pqlctsizeFromFormatAndSize;
        
end