function img = filterEntropy(img, ksize)
%
% img = filterEntropy(img, ksize)
%
% description:
%    local entropy filter
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

%img  = img - mean(img(:));
img = entropyfilt(img, ksize);

end