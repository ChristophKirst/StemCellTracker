function [shift, quality] = align2ImagesOnGridByHugin(imgs, varargin)
%
% shift = align2ImagesOnGridByHugin(imgs, param)
%
% description:
%     aligns two images using overlap weighted cross correlation calculated via fft
%     assuming images are prealigned as indicated by the input cell {left; right} {bottom, up}, {{down}, {top}}
%
% input:
%     imgs         images to align as prealigned cell array
%     param        parameter struct with entries
%                  as in align2ImagesLeftRightByCorrelation
%                  as in align2ImagesLeftRightParameter
%
% output:
%     shift        the shift between origin the images imgs in pixel coordinates and pixel units
%     quality      (optional) quality measure of alignment
%
% See also: align2ImagesLeftRightParameter, align2ImagesLeftRightByHugin

[img1, img2, per] = align2ImagesLeftRightOrient(imgs);

[shift, quality] = align2ImagesLeftRightByHugin(img1, img2, varargin{:});

shift = shift(per);

end
         


