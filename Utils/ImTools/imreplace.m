function repl = imreplace(image, subimage, coords, chop)
%
% repl = imreplace(image, subimage, coords)
%
% description: 
%    replaces a subimage in image at coordinates coords
%
% input: 
%    image     original image
%    subimage  replacement sub image
%    coords    p,q,l corner coordinates from which to start replacement
%    chop      if subimage does not fit then chop if true, otherwise gerneate error
%
% output:
%    repl      image in which subimage is replaced
%
% See also: imextract, imfind

if nargin < 4
   chop = false;
end
if ischar(chop)
   if strcmp(chop, 'chop')
      chop = true;
   else 
      chop = false;
   end
end

isize = size(image);
ssize = size(subimage);

if numel(coords) ~= ndims(image) || numel(coords) ~= ndims(subimage) || any(coords > isize)
   disp(coords)
   ndims(image)
   ndims(subimage)
   error('imreplace: inconsistent image dimensions or positions!')
end

if any(ssize + coords - 1 > isize)
   if chop
      ssize = min(ssize, isize - coords + 1);
      si = cell(1,length(ssize));
      for i = 1:length(ssize)
         si{i} = 1:ssize(i);
      end
      subimage = subimage(si{:});
   else
      error('imreplace: subimage to large!')
   end
end

repl = image;
if isempty(subimage)
   return
end

idx = arrayfun(@(i,j)(i:j), coords , coords + ssize-1, 'UniformOutput', false);
repl(idx{:}) = subimage;

end




