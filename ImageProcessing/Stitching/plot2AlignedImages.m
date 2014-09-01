function plot2AlignedImages(img1, img2, shift)
%
% plot2AlignedImages(img1, img2, shift)
%

if ~iscell(shift)
   shift = {zeros(1,ndims(img1)), shift};
end

plotAlignedImages({img1, img2}, shift)

end