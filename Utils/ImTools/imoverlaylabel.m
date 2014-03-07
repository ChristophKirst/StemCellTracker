function iol = imoverlaylabel(img, label, varargin)
%
% iol = imoverlaylabel(img, label, param)
%
% descriptiopn:
%    takes grayscale image and overlays colorized label
%
% input:
%    img      grayscale image
%    label    labeled image
%    param    (optional) paramter struct for color map as in imlabelcolormap
%
% output:
%    iol       overlay of grayscale with colorized label
%
% See also: imcolorize, imoverlay, imlabelcolormap

iol = gray2rgb(img);
iol = iol / max(iol(:));
iol = reshape(iol, [],3);
imgcl = imcolorize(label, varargin{:});
imgcl = reshape(imgcl, [],3);
idx = find(label);
iol(idx,:) = imgcl(idx,:);
iol = reshape(iol, [size(img) 3]);

end