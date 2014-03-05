function box = imextractbox(image, center, ksize, padding)
%
% ie = imextractbox(image, pos, ksize, padding)
%
% description: 
%     extracts the region specified by rect in the image
%
% input:
%     image    the image
%     pos      center pixel coordinates
%     ksize    box size (h x w (x l)
%     padding  (optional) '
%
% output:
%     box       sub image
%
% See also: imcrop, imextract

if nargin < 4
   padding = 'none';
end

[ol, or] = filteroffsets(ksize);

sub2ind

rectlo = center - ol;
recthi = center + or;

isize = size(image);

switch padding
   case 'none'
      rectlo(rectlo < 1) = 1;
      idx = recthi > isize;
      recthi(idx) = isize(idx);

   otherwise
      
      if any(rectlo < 1) || any(recthi > isize)
         offset = max(ol, or);
         image = padarray(image, offset, padding);
         rectlo =rectlo + offset; recthi = recthi + offset;
      end
end

for i = ndims(image):-1:1
   rect{i} = rectlo(i):recthi(i);
end

box = image(rect{:});

end
