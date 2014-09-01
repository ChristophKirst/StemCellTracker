function [shift, quality] = align2ImagesLeftRightByHugin(img1, img2, varargin)
%
% [shift, quality] = align2ImagesLeftRightByHugin(img1, img2, param)
%
% description: 
%     aligns two images using hugin feature key alignment
%     assuming images are prealigned (img1 - left, img2 - right)
%
% input:
%     img1,img2    images to align
%     param        parameter struct with entries
%                  as in align2ImagesLeftRightParameter
%                  as in parameter as in hialign
%
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units
%     quality      (optional) quality measure of alignment = 1
%
% See also: align2ImagesLeftRightParameter, hialign

param = parseParameter(varargin{:});
[~, s1, ~, minovl, maxovl, maxshift] = align2ImagesLeftRightParameter(img1, img2, param);

% cut relavant regions
img1c = extract(img1, s1(1)-maxovl:s1(1), 1);
img2c = extract(img2, 1:maxovl, 1);

%figure(96); clf; implottiling({img1c, img2c});

% align via hugin

shift = align2ImagesByHugin(img1c, img2c, varargin{:});

% correct for image cutting and out of bounds shifts
ms = [maxovl, maxshift];
id = shift > ms;
shift(id) = ms(id);

ms = [minovl, -maxshift];
id = shift < ms;
shift(id) = ms(id);

shift(1) = shift(1) + s1(1) - maxovl - 1;
quality = 1;

end