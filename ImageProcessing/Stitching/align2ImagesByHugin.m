function [shift, quality] = align2ImagesByHugin(img1, img2, varargin)
%
% [shift, quality] = align2ImagesByHugin(img1, img2, param)
%
% description:
%     globally aligns two images img1, img2 using Hugin panorama tool
%
% input:
%     img1, img1   images to be alligned, assumed to have same orientation
%     param        (optional)  struct with entries as in hialign function
%
% output:
%     shift        the shift between origin of img1 to orgigin img2 in pixel coordinates and pixel units
%     quality      (optional) quality measure of the alignment = 1
%
% See also:  hialign

shift = alignImagesByHugin({img1, img2}, varargin{:});
shift = shift{2};
quality = 1;

end
         


