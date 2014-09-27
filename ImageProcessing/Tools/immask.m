function masked = immask(image, mask)
%
% masked = immask(image, mask)
%
% description:
%      applies mask to image
%
% input:
%     image    image to mask
%     mask     the maks
%
% output:
%    masked    masked image

masked = image .* cast(mask, 'like', image);

end