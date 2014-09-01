function [shift, quality] =  align2ImagesOnGridByOptimization(imgs, varargin)
%
% img = align2ImagesByOptimization(imgs, param)
%
% description:
%     aligns two images using optimization
%     assuming images are prealigned as indicated by the input cell {left; right} {top, bottom}, {{down}, {up}}
%
% input:
%     imgs         images to align as prealigned cell array
%     param        parameter struct with entries
%                  as in align2ImagesLeftRightByOptimization
%                  as in align2ImagesLeftRightParameter
%
% output:
%     shift        the shift between of the images imgs in pixel coordinates and pixel units
%     quality      (optional) quality measure of alignment
%
% See also: align2ImagesLeftRightParameter, align2ImagesLeftRightByOptimization

[img1, img2, per] = align2ImagesLeftRightOrient(imgs);

[shift, quality] =  align2ImagesLeftRightByOptimization(img1, img2, varargin{:});

shift = shift(per);

end