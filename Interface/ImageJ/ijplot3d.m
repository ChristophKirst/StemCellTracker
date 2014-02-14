function universe = ijplot3d(image, varargin)
%
% universe = ijplot3d(image, varargin)
%
% description:
%    starts ImageJ interface
%
% input:
%    image      image to plot
%    varargin   parameter 'ParameterName' ParameterValue
%               'ImageSequence'   (true = display image as sequence)
%               'PixelDepth'      (z depth of pixel w.r.t to x,y)
%
% output:
%    universe    Universe object representing the 3D viewer window
%
% See also: ijstart


param.sequence = false;
param.pixel_depth = 1;

for n = 1:2:length(varargin)
   switch(lower(varargin{n}))
      case 'imagesequence'
         param.sequence = varargin{n+1};
      case 'pixeldepth'
         param.pixel_depth = varargin{n+1};
   end
end

% convert image to [X Y Z C]
% matlab uses [X Y C Z]
dim = ndims(image);

if dim == 4 % color image
   image = permute(image, [1 2 4 3]);
else % stack of gray scale images
   image = cat(4, image,image,image);
end
% need to convert to uint8 ?

% create image

% get ij instance



imp = MIJ.createColor('Figure', image, param.sequence);

% pixel depth
calibration = ij.measure.Calibration();
calibration.pixelDepth = param.pixel_depth;
imp.setCalibration(calibration);

%new 3D viewer
universe = ij3d.Image3DUniverse();
universe.show();
universe.addVoltex(imp);

end