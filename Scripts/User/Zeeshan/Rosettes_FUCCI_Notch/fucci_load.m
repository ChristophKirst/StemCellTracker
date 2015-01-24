function img = fucci_load(is, time, pos)

ds = is.dataSize;
img = zeros(ds(1),ds(2),3);

% Scale
scalemax = [4000  360     3400];
scalemin = [100   130     150];

% load image
imgB = is.data('T', time, 'S', pos, 'C', '4CFP');
imgR = is.data('T', time, 'S', pos, 'C', '3RFP');
imgG = is.data('T', time, 'S', pos, 'C', '2YFP');

img(:,:,1) = imgR; img(:,:,2) = imgG; img(:,:,3) = imgB;

for c = 1:3
   img(:,:,c) = (img(:,:,c) - scalemin(c))/(scalemax(c) - scalemin(c));
   %img(:,:,c) = mat2gray(img(:,:,c));
end

% enhance red
for c = 1
   %img(:,:,c) = 1.5 * img(:,:,c);
end


end