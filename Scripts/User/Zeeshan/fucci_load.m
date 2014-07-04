function img = fucci_load(exp, time)

channelnames =   {'w2CFP',  'w3RFP',  'w4YFP', 'w1DIC'};

img = zeros( 1344,1024,3);

% Scale
scalemax = [4000  360     3400];
scalemin = [100   130     150];

% load image
imgB = imread(exp.ImageFile('time', time, 'channel', channelnames{1}))';
imgR = imread(exp.ImageFile('time', time, 'channel', channelnames{2}))';
imgG = imread(exp.ImageFile('time', time, 'channel', channelnames{3}))';

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