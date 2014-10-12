function img = filterRange(img, ksize)
%
% img = filterRange(img, ksize)
%
% description:
%    local range filter
%
% input:
%    img          image to be filtered
%    ksize        h x w (xl) kernel size, ones(ksize) or mask
%
% output:
%    img          filtered image

if length(ksize) == numel(ksize)
   ksize = ones(ksize);
end

img = rangefilt(img, ksize);

end