function c = imarisrgb2spectrum(rgb)
%
% c = imarisrgb2spectrum(rgb)
%
% converst a rgb color to best fit 



%%
vRed = 1.00;
vGreen = 1.0;
vBlue = 0.6;
vAlpha = 0;
%disp(sprintf('Range 0-1: Red = %f, green = %f, blue = %f, alpha = %f', vRed, vGreen, vBlue, vAlpha))
vRGBA = [vRed, vGreen, vBlue, vAlpha];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]) % combine different components (four bytes) into one integer

%%
sf.SetColorRGBA(vRGBA)

%%
sf.SetSelectedIndices([])

%%
sf.GetSelectedIndices()

