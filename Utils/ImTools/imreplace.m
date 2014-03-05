function repl = imreplace(image, subimage, coords)
%
% repl = imreplace(image, subimage, coords)
%
% description: 
%    replaces a subimage in image at coordinates coords
%
% input: 
%    image     original image
%    subimage  replacement sub image
%    coords    h,w,l corner coordinates from which to start replacement
%
% output:
%    repl      image in which subimage is replaced
%
% See also: imextract, imfind

isize = size(image);
ssize = size(subimage);

if numel(coords) ~= ndims(image) || numel(coords) ~= ndims(subimage)
   coords
   ndims(image)
   ndims(subimage)
   error('imreplace: inconsistent image dimensions!')
end

if any(ssize + coords - 1 > isize)
   error('imreplace: subimage to large!')
end

repl = image;
idx = arrayfun(@(i,j)(i:j), coords , coords + ssize-1, 'UniformOutput', false);
repl(idx{:}) = subimage;

end




