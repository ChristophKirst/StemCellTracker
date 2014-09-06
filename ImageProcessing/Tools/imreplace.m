function img = imreplace(img, subimg, coords, chop)
%
% repl = imreplace(img, subimg, coords)
%
% description: 
%    replaces a subimage in image at coordinates coords
%
% input: 
%    img     original image
%    subimg  replacement sub image
%    coords  p,q,l corner coordinates from which to start replacement
%    chop    if subimage does not fit then chop if true, otherwise gerneate error
%
% output:
%    repl      image in which subimage is replaced
%
% See also: imextract, imfind

if nargin < 4
   chop = false;
else
   if ischar(chop) && strcmp(chop, 'chop')
      chop = true;
   else
      chop = false;
   end
end

isize = size(img);
ssize = size(subimg);
dim =  ndims(img);

if numel(coords) ~= dim || dim ~= ndims(subimg) || any(coords > isize)
   error('imreplace: inconsistent image dimensions or positions!')
end

if any(ssize + coords - 1 > isize)
   if chop
      ssize = min(ssize, isize - coords + 1);
      si = cell(1,length(ssize));
      for i = 1:dim
         si{i} = 1:ssize(i);
      end
      subimg = subimg(si{:});
   else
      error('imreplace: subimage to large!')
   end
end

if isempty(subimg)
   return
end

idx = arrayfun(@(i,j)(i:j), coords , coords + ssize-1, 'UniformOutput', false);
img(idx{:}) = subimg;

end




