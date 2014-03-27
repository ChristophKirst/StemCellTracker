function c = imarisrgb2spectrum(r,g,b)
%
% c = imarisrgb2spectrum(rgb)
%
% description:
%     convert a rgb color to closest color on spectrum (inverse Hue) color
%
% input:
%     rgb     rgb values [r,g,b]
%
% output:
%    c        color indices for use with 'Spectrum' color map in Imaris
%
% See also: imarisrgb2color, imarisrgb2statistics

if nargin > 1
   rgb = [r,g,b];
else
   rgb = r;
end

%rgb = mat2gray(rgb);
c = rgb2hsv(rgb);

cmax = 14./16.;
c(c>cmax) = cmax;
c = (cmax - c(:,1))/cmax;

end

%color coding in imaris
% vRed = 1.00;
% vGreen = 1.0;
% vBlue = 0.6;
% vAlpha = 0;
% %disp(sprintf('Range 0-1: Red = %f, green = %f, blue = %f, alpha = %f', vRed, vGreen, vBlue, vAlpha))
% vRGBA = [vRed, vGreen, vBlue, vAlpha];
% vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
% vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]) % combine different components (four bytes) into one integer
% obj.SetColorRGBA(vRGBA)

