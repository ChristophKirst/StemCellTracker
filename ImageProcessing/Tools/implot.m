function varargout = implot(img, varargin)
%
% implot(img, varargin)
%
% description:
%    implot plots 2d or 3d image
%
% input:
%    img     image to plot
%
% See also: implot2d, implot3d

% switch 2d / 3d image

imf = lower(imfrmtFormat(img));

switch imf
   case {'xy', 'xyc'}
      img = implot2d(img, varargin{:});
   case {'xyz', 'xyzc'}
      img = implot3d(img, varargin{:});
   otherwise
      error('implot: image format is not xy, xyc, xyz, xyzc but %s', imf);
end


if nargout > 0
   varargout{1} = img;
end
     
