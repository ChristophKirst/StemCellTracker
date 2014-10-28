function img = imreplace(img, subimg, coords, varargin)
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

if isempty(subimg)
   return
end

param = parseParameter(varargin);
chop = getParameter(param, 'chop', false);

isize = size(img);
ssize = size(subimg);
dim =  ndims(img);

if numel(coords) ~= dim || dim ~= ndims(subimg)
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
      ssize = size(subimg);
      if isempty(subimg)
         return
      end 
   else
      error('imreplace: subimage to large!')
   end
end

if any(coords < 1)
   if chop
      si = cell(1,length(ssize));
      istart = max(2-coords, 1);
      for i = 1:dim
         si{i} = istart(i):ssize(i);
      end
      subimg = subimg(si{:});
      coords = max(1, coords);
      ssize = size(subimg);
      
      if isempty(subimg)
         return
      end 
      
   else
      error('imreplace: subimage coordinates out of range!')
   end
end


if isempty(subimg)
   return
end

idx = arrayfun(@(i,j)(i:j), coords , coords + ssize-1, 'UniformOutput', false);
img(idx{:}) = subimg;

end




