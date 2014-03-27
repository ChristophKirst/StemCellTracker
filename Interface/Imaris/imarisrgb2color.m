function c = imarisrgb2color(r, g, b, a)
%
% c = imarisrgb2color(rgb)
%
% description:
%     convert a rgb color to color integer used by imaris
%
% input:
%     rgb     rgb values [r,g,b] between 0,1
%
% output:
%    c        color integer
%
% See also: imarisrgb2specturm

if nargin == 1
   g = r(:,2);
   b = r(:,3);
   r = r(:,1);
   
   if size(rgb,2) == 4
      a = r(:,4);
   else
      a = zeros(size(rgb,1),1);
   end
   
elseif nargin == 3
   a = zeros(size(r,1),1);
elseif nargin == 2
   error('imarisrgb2color: expects r,g,b, (a) or rgb(a) input');
end

% %disp(sprintf('Range 0-1: Red = %f, green = %f, blue = %f, alpha = %f', vRed, vGreen, vBlue, vAlpha))

c = [r, g, b, a];
c = round(c * 255); % need integer values scaled to range 0-255
c = uint32(c * [1; 256; 256*256; 256*256*256]); % combine different components (four bytes) into one integer

% obj.SetColorRGBA(vRGBA)

end

