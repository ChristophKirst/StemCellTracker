function [ie, rectlow, recthigh] = imextract(image, rectlow, recthigh)
%
% ie = imextract(image, rect)
% ie = imextract(image, rectlow, recthigh)
%
% description: 
%     extracts the region specified by rect in the image
%
% input:
%     image    the image
%     rect     rectangle in pixel coordinates [rectlow recthigh] or 'BoundingBox'
%     rectlow  lower corner
%     recthigh higher corner
%
% output:
%     ie       cropped image
%
% See also: imcrop

if isequal(rectlow, 'BoundingBox')
   [rectlow, recthigh] = imboundingbox(image);
   d = length(rectlow);
elseif nargin < 3
   d = length(rectlow) / 2;
   recthigh = rectlow(d+1:2*d);
   rectlow = rectlow(1:d);
else 
   d = length(rectlow);
end

if d ==2
   ie = image(rectlow(1):recthigh(1), rectlow(2):recthigh(2));
else
   ie = image(rectlow(1):recthigh(1), rectlow(2):recthigh(2), rectlow(3):recthigh(3));
end
   
   
