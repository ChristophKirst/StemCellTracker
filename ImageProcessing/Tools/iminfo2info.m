function i = iminfo2info(imfi)
%
% i = iminfo2info(d)
%
% description:
%    infers image information from the struct returned by imfinfo and converts it to a ImageInfo class
%    it gives consistent data when readin gthe image with imread
%
% input:  
%   imfi   struct returned by imfinfo
%   i      ImageInfo class
%

i = ImageInfo();
i.idatasize   = [imfi.Height, imfi.Width];
i.idataformat = imsize2format(i.isize);
i.idataclass  = ['unit', imfi.BitsPerSample];
i.irawformat  = i.idataformat;
i.irawsize    = i.idatasize;

i.pqlctsizeFromFormatAndSize;

switch imfi.ColorType
   case 'grayscale'
      i.icolor = {'gray'};
   case 'truecolor'
      i.icolor = {'r', 'g', 'b'};
   case 'indexed'
      cl = 1;
      id = find(i.iformat == 'c', 1, 'first');
      if ~isempty(id)
         cl = i.isize(id);
      end
      i.icolor  = imcolorlist(cl);
end
        
end