function swap = imswapcoords(image)
% 
% swap = imswapcoords(image)
%
% description: 
%      swaps between pixel and space coordinates
%
% input:
%      image   image to be swapped
%
% output:
%      swap     imaged with swapped coordinates 
%
% See also: imformat, imfind3d


dim = ndims(image);
if dim < 2
   error('imswapcoords: expects image as input');
end

ids = 1:ndims(image);
ids(1:2) = [2,1];

swap = image(ids);

end

