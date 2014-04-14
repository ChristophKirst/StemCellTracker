function impl = implus(image, subimage, coords, chop)
%
% repl = implus(image, subimage, coords, chop)
%
% description: 
%    adds the subimage to image at coordinates coords
%
% input: 
%    image     original image
%    subimage  subimage to be added at certain location
%    coords    p,q,l corner coordinates in image of subsimage
%    chop      if subimage does not fit then chop if true, otherwise gerneate error
%
% output:
%    impl      image in which subimage is added at position coodrs 
%
% See also: imreplace, imextract, imfind

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
   disp coords:
   disp(coords)
   disp ndims(image):
   ndims(image)
   disp ndims(subimage):
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

impl = image;
if isempty(subimage)
   return
end

idx = arrayfun(@(i,j)(i:j), coords , coords + ssize-1, 'UniformOutput', false);
impl(idx{:}) = impl(idx{:}) + subimage;

end




