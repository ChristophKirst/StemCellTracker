function iol = imoverlaylabel(img, label)
%
% iol = imoverlaylabel(img, label)
%
% descriptiopn:
%    takes grayscale image and overlays colorized label
%
% input:
%    img   grayscale image
%    label labeled image
%
% output:
%    iol   overlay of grayscale with colorized label
%
% See also: imcolorize, imoverlay

iol = gray2rgb(img);
iol = reshape(iol, [],3);
imgcl = imcolorize(label);
imgcl = reshape(imgcl, [],3);
idx = find(label);
iol(idx,:) = imgcl(idx,:);
iol = reshape(iol, [size(img) 3]);

end