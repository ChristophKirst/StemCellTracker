function iinfo = imfinfo2info(iminfo)
%
%  info = imfinfo2info(iinfo)
%
% description:
%   convert from info obtained form matlabs imfinfo to ImageInfo class


iinfo = ImageInfo();

iinfo.idatasizePQLCT = [iminfo.Width, iminfo.Height, 1, iminfo.SamplesPerPixel, 1];

if iminfo.SamplesPerPixel == 1
   iinfo.isize      = [iminfo.Width, iminfo.Height];
   iinfo.iformat    = ['pq'];
else
   iinfo.isize      = [iminfo.Width, iminfo.Height, iminfo.SamplesPerPixel];
   iinfo.iformat    = ['pqc'];
end

switch iminfo.ColorType
   case 'truecolor'
      iinfo.icolor = imcolorlist(3);
   case 'grayscale'
      iinfo.icolor = {'gray'};
   otherwise
      iinfo.icolor = imcolorlist(iinfo.sizeC);
end

bps = max(iminfo.BitsPerSample);
switch bps
   case 8
      iinfo.iclass = 'uint8';
   case 16
      iinfo.iclass = 'uint16';
   case 32
      iinfo.iclass = 'uint32';   
   case 64
      iinfo.iclass = 'uint64';
   otherwise
      iinfo.iclass = 'double';
end

iinfo = determineScale(iinfo, iminfo);

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


   %TODO:  %ZVI  %LIF
end

