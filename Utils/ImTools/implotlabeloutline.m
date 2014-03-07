function implotlabeloutline(image, label, varargin)
%
% implotlabeloutline(image, label, param)
%
% description:
%     plots image as tiling together with pixel surface of the labeled image
%
% input:
%     image    grayscale image
%     label    labeled image
%     param    (optional) paramter struct for color map as in imlabelcolormap
%
% See also: impixelsurface, imsurfaceplot3d, imoverlay, imlabelcolormap

implottiling(imoverlaylabel(image, impixelsurface(label), varargin{:}))

end