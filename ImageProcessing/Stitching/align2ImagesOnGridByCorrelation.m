function [shift, quality] = align2ImagesOnGridByCorrelation(imgs, varargin)
%
% [shift, quality] = align2ImagesOnGridByCorrelation(imgs, param)
%
% description: 
%     aligns two images using overlap weighted cross correlation calculated via fft
%     assuming images are prealigned as indicated by the input cell {left; right} {top, bottom}, {{down}, {up}}
%
% input:
%     imgs         images to align as cell array
%     param        parameter struct with entries
%                  as in align2ImagesLeftRightByCorrelation
%                  as in align2ImagesLeftRightParameter
%
% output:
%    shift        the shift between origin of the images in pixel coordinates and pixel units
%    quality      (optional) quality measure of alignment
%
% See also: align2ImagesLeftRightParameter, align2ImagesLeftRightByCorrelation

[img1, img2, per] = align2ImagesLeftRightOrient(imgs);

[shift, quality] = align2ImagesLeftRightByCorrelation(img1, img2, varargin{:});

shift = shift(per);

end