function img = immask(img, imgmask)
%
% masked = immask(image, mask)
%
% description:
%      applies mask to image
%
% input:
%     img      image to mask
%     mask     the maks
%
% output:
%     img    masked image

if ndims(img) == ndims(imgmask)
   img = img .* cast(imgmask, 'like', img);
elseif ndims(img) == ndims(imgmask)+1
   si = size(img);
   si(1:length(si)-1) = 1;
   si = num2cell(si);
   imgmask = repmat(imgmask, si{:});
   img = img .* cast(imgmask, 'like', img);
end