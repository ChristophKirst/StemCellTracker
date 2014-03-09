function [ie, rectlow, recthigh] = imextract(image, rectlow, recthigh)
%
% ie = imextract(image, rect)
% ie = imextract(image, rectlow, recthigh)
%
% description: 
%     extracts the region specified by rect in the image using pixel coordinates
%
% input:
%     image    the image
%     rect     rectangle in pixel coordinates [rectlow recthigh] or 'BoundingBox'
%     rectlow  lower corner (in p,q,(l) coordinates)
%     recthigh higher corner (in p,q,(l) coordinates)
%
% output:
%     ie       extracted subimage 
%
% See also: imcrop, imextractbox

d = ndims(image);

if isequal(rectlow, 'BoundingBox')
   [rectlow, recthigh] = imboundingbox(image);
elseif nargin < 3
   recthigh = rectlow(d+1:2*d);
   rectlow = rectlow(1:d);
end

for i = d:-1:1
   rect{i} = rectlow(i):recthigh(i);
end

ie = image(rect{:});

end

   
   
