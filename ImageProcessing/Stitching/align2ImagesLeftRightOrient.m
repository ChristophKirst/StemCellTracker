function [img1, img2, per] = align2ImagesLeftRightOrient(imgs)
%
% [img1, img2, per] =  align2ImagesLeftRightOrient(imgs)
%
% descritpion:
%   orient two images for global image aligment by transformaing into left right aligment
%
% input:
%   imgs       images to align
%           
% output:
%   img1,img2   oriented images
%   per         permutation to get back to original image

if numel(imgs) ~= 2
   error('align2ImagesLeftRightOrient: expect precisle cell array with two images!');
end
img1 = imgs{1}; img2 = imgs{2};

% orientation

si = size(imgs);
dim = ndims(imgs);
pos = find(si == 2, 1);
per = 1:dim;
per(1) = pos; per(pos) = 1;
img1 = permute(img1, per);
img2 = permute(img2, per);

end







