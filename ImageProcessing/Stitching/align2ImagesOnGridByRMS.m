function [shift, quality] = align2ImagesOnGridByRMS(imgs, varargin)
%
% img = align2ImagesOnGridByRMS(imgs, param)
%
% description:
%     aligns two images using overlap weighted root means squares calculated via fft
%     assuming images are prealigned as indicated by the input cell {left; right} {top, bottom}, {{down}, {up}}
%
% input:
%     imgs         images to align as prealigned cell array
%     param        parameter struct with entries
%                  as in align2ImagesLeftRightByRMS
%                  as in align2ImagesLeftRightParameter
%
% output:
%     shift        the shift between origin of the images imgs in pixel coordinates and pixel units
%     quality      (optional) quality measure of alignment
%
% See also: align2ImagesLeftRightParameter, align2ImagesLeftRightByRMS

[img1, img2, per] = align2ImagesLeftRightOrient(imgs);

[shift, quality] = align2ImagesLeftRightByRMS(img1, img2, varargin{:});

shift = shift(per);

end