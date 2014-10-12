function img = filterStd(img, ksize)
%
% img = filterStd(img, ksize)
%
% description:
%    local standard deviation filter
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

img = stdfilt(img, ksize);

end