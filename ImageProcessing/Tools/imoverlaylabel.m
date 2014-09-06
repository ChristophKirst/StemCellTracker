function iol = imoverlaylabel(img, label, intensity, varargin)
%
% iol = imoverlaylabel(img, label, intensity, param)
%
% descriptiopn:
%    takes grayscale image and overlays colorized label
%
% input:
%    img       grayscale image
%    label     labeled image
%    intensity (optional) keep intensity of img in the colorized label (true)
%    param     (optional) paramter struct for color map as in imlabelcolormap
%
% output:
%    iol       overlay of grayscale with colorized label
%
% See also: imcolorize, imoverlay, imlabelcolormap

if nargin < 3 || isempty(intensity)
   intensity = true;
end


iol = gray2rgb(img);
iol = iol / max(iol(:));
iol = reshape(iol, [],3);
imgcl = imcolorize(label, varargin{:});
imgcl = reshape(imgcl, [],3);
idx = find(label);

if intensity
   iol(idx,:) = imgcl(idx,:) .* iol(idx,:);
else
   iol(idx,:) = imgcl(idx,:);
end

iol = reshape(iol, [size(img) 3]);

end