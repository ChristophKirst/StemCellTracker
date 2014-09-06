function values = imfiltervalues(image, mask, filter)
%
% values = imfiltervalues(image, mask, filter)
%
% description: 
%     calculates for the white pixels in mask the value using the
%     fitler. mask can also be a 1xn or 2xn matrix of indices or pixel coordinates
%
% input:
%     image    the image
%     mask     binariy image indicating the pixels to aclulate means for
%              or list of pixel indices [i1, i2, ...]
%              or list of pixel coordinates [p1, p2, ...; q1, q2, ...;(l1, l2 ,...)]        
%     filter   the kernel matrix for the filter
%
% output:
%     values   list of filter values the pixels in mask
%
% note:
%     on the border the fitler is simply cut off
%     if coordinates are given these are pixel coordinates
%
% See also: roifilt2


if nargin < 3
   filter = 3;
end

d = ndims(image);
if d ~=2 && d ~= 3
   error('imfiltervalues: image dimensions must be 2 or 3!')
end

if size(mask,1) == 1
   ind = mask;
   mask = zeros(size(image));
   mask(ind) = 1;
elseif size(mask,1) == 2 && d == 2
   ind = mask;
   mask = zeros(size(image));
   mask(sub2ind(size(image), ind(1,:), ind(2,:))) = 1;
elseif size(mask,1) == 3 && d == 3
   ind = mask;
   mask = zeros(size(image));
   mask(sub2ind(size(image), ind(1,:), ind(2,:), ind(3,:))) = 1;
end

if ~isequal(size(mask), size(image))
   error('imfiltervalues: mask and image dimensions disagree!')
end

if length(filter) == 1
   if filter == 0
      values = image(mask > 0);
      return
   end

   if d == 2
      filter = fspecial2('disk', filter);
      %filter = filter / sum(filter(:));
   else
      filter = fspecial3('disk', filter);
   end
end


switch d
   case 2

      values = roifilt2(filter, image, mask);
      values = values(mask > 0);
      
   case 3
      
      pos = find(mask);
      
      if length(pos) > numel(image)/10 %if lots of points fft is faster
         image = imfilter(image, filter, 'symmetric');
         values = image(pos);
      else
         
         [oh, ow, ol] = size(filter);
         [ohl, ohr] = filteroffsets(oh);
         [owl, owr] = filteroffsets(ow);        
         [oll, olr] = filteroffsets(ol);

         [ih,iw,il] = ind2sub(size(image), pos);
         
         image = padarray(image, [oh, ow, ol] , 'symmetric');

         for i = length(pos):-1:1
            values(i) = sum(sum(sum( ...
                   filter .* image((oh + ih(i) - ohl) : (oh + ih(i) + ohr),...
                                   (ow + iw(i) - owl) : (ow + iw(i) + owr),...
                                   (ol + il(i) - oll) : (ol + il(i) + olr))...
                                   )));
         end
      end
end

end




