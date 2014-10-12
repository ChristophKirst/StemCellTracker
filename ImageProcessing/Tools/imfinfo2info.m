function iinfo = imfinfo2info(finfo)
%
%  info = imfinfo2info(finfo)
%
% description:
%   convert info obtained form imfinfo to ImageInfo class
%
% input:
%   finfo    info struct obtained by imfinfo

iinfo = ImageInfo();

if finfo.SamplesPerPixel == 1
   iinfo.idatasize      = [finfo.Width, finfo.Height];
   iinfo.idataformat    = 'pq';
else
   iinfo.idatasize      = [finfo.Width, finfo.Height, finfo.SamplesPerPixel];
   iinfo.idataformat    = 'pqc';
end

iinfo.pqlctsizeFromFormatAndSize;

switch finfo.ColorType
   case 'truecolor'
      iinfo.icolor = imcolorlist(3);
   case 'grayscale'
      iinfo.icolor = {'gray'};
   case 'indexed'
      iinfo.icolor  = imcolorlist(iinfo.datasizeC);
   otherwise
      iinfo.icolor = imcolorlist(iinfo.datasizeC);
end

bps = max(finfo.BitsPerSample);
switch bps
   case 8
      iinfo.idataclass = 'uint8';
   case 16
      iinfo.idataclass = 'uint16';
   case 32
      iinfo.idataclass = 'uint32';   
   case 64
      iinfo.idataclass = 'uint64';
   otherwise
      iinfo.idataclass = 'double';
end

iinfo = determineScale(iinfo, finfo);

end



%infer spatial scales from meta data, sets: scaleP, scaleQ, scaleL, times
function iinfo = determineScale(iinfo, iminfo)
   % extract info from meta data:

   % TIF XYResolution is number of Pixel per ResolutionUnit
   if strcmp(iminfo.Format, 'tif')
      xr = iminfo.XResolution;
      yr = iminfo.YResolution;

      iinfo.iscale = [yr, xr];
      %iinfo.iscaleFormat = 'pq';
      iinfo.iunit =  iminfo.ResolutionUnit;

      return
   end
   
   %Todo: bmp, png etc
end

