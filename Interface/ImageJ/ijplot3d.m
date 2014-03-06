function universe = ijplot3d(image, varargin)
%
% universe = ijplot3d(image, varargin)
%
% description:
%    plots a 3d image using the ImageJ 3DViewer
%
% input:
%    image      3D image to plot
%    varargin   parameter 'ParameterName' ParameterValue
%               'PixelDepth'      z depth of pixel w.r.t to x,y (1)
%               'Name'            name of image ('Figure')
%               'Universe'        add the image to an exsiting universe
%
% output:
%    universe    Universe object representing the 3D viewer window
%
% See also: ijstart

pixel_depth = 1;
name = '3DImage (Matlab)';
univ = [];

for n = 1:2:length(varargin)
   switch(lower(varargin{n}))
      case 'pixeldepth'
         pixel_depth = varargin{n+1};
      case 'name'
         name = varargin{n+1};
      case 'universe'
         univ = varargin{n+1};   
   end
end


dim = ndims(image);

if dim == 4 % color image
   % nothing to be done / maybe check format hwlc !!
else % stack of gray scale images
   image = cat(4, image,image,image);
end
% convert to uint8 

image = cast(255 * image / max(image(:)), 'uint8');

imp = MImageJ.createImage(name, 'hwlc', image)

% pixel depth
calibration = ij.measure.Calibration();
calibration.pixelDepth = pixel_depth;
imp.setCalibration(calibration);

% plot
if isa(univ, 'ij3d.Image3DUniverse')
   universe = univ;
else
   universe = ij3d.Image3DUniverse();
end
universe.show();
universe.addVoltex(imp);

end