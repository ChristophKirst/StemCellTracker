function profile = imray(image, startp, endp, n)
%
% profile = imray(image, startp, endp, n)
%
% description:
%    returns n profile points along the ray from start point to end point
%    using pixel coorindates [p,q,l]
%
% input:
%    image     2d or 3d gray scale image
%    startp    starting point in pixel coordinates [p; q(; l), ...]
%    endp      end point in pixel coordinates  [p; q(; l), ...]
%    n         number of points
%
% output:
%    profile   n intensity points along ray [i1; i2; ...]
%
% See also: interp3, improfile

d = ndims(image);
if size(startp,1) ~= d || size(endp,1) ~= d 
   error('imray: dimension of image and ray coordinates does not match!')
end
k = size(startp,2);
if size(endp,2) ~= k 
   error('imray: numer of ray start and end points does not match!')
end 

if nargin < 4
   n = max(round(sqrt(sum((startp-endp).^2,1))));
end
   
profile = zeros(n,k);


for i=k:-1:1
   switch d
      case 2
         profile(:,i) = improfile(image, [startp(2,i) endp(2,i)], [startp(1,i) endp(1,i)], n); % [p,q] = [y,x], [x,y] = [q,p] in matlab improfile !!
      
      case 3
         %create interpolation points
         xp = linspace(startp(2,i), endp(2,i), n);
         yp = linspace(startp(1,i), endp(1,i), n);
         zp = linspace(startp(3,i), endp(3,i), n);

         profile(:,i) = interp3(double(image), xp, yp, zp);      
   end
end

end



