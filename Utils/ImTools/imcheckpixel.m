function [p, id] = imcheckpixel(p, isize)
%
% p = imcheckpixel(p, isize)
%
% description:
%   checks if pixel coordinates p are in the image
%
% input:
%   p      pixel coordinates [x,y(,z)] 
%   isize  image size
%
% output:
%   p      pixel that are in the image

dim = size(p,2);
id = ones(size(p,1),1);

isize = padright(isize, dim, 1);

for d = 1:dim
   id = and(id, 1 <= p(:,d));
   id = and(id, p(:,d) <= isize(dim));
end
p = p(id, :);

end