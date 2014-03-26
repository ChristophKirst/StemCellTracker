function c = imarisrgb2spectrum(rgb)
%
% c = imarisrgb2spectrum(rgb)
%
% description:
%     convert a rgb color to closest color on spectrum (ínverse Hue) color
%
% input:
%     rgb     rgb values [r,g,b]
%     ncols   (optional) number of different colors (256)
%
% output:
%    c        color indices fur use with spectrum map
%
% See also: imarisrgb2color, hsv, rgb2hsv

rgb = mat2gray(rgb);
c = rgb2hsv(rgb);
c = 1 - c(:,1);

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

