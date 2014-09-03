function [shift, quality] = align2ImagesLeftRightByCorrelation(img1, img2, varargin)
%
% [shift, quality] = align2ImagesLeftRightByCorrelation(img1, img2, param)
%
% description: 
%     aligns two images using overlap weighted cross correlation calculated via fft
%     assuming images are prealigned (img1 - left, img2 - right)
%
% input:
%     img1,img2    images to align
%     param        parameter struct with entries
%                  .mean    if true subtract mean (true)
%                  as in align2ImagesLeftRightParameter
%
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units
%
% See also: align2ImagesLeftRightParameter

param = parseParameter(varargin{:});
[dim, s1, s2, minovl, maxovl, maxshift] = align2ImagesLeftRightParameter(img1, img2, param);

% cut relavant regions
img1c = extract(img1, s1(1)-maxovl:s1(1), 1);
img2c = extract(img2, 1:maxovl, 1);

%pad images if it is not the same size
s1c = size(img1c); s2c = size(img2c);
sic = max(s1c,s2c);
img1c = padarray(img1c, sic - s1c, 0, 'post');
img2c = padarray(img2c, sic - s2c, 0, 'post');


% pad border to ensure no overlaps for maximal shifts
pm = [sic(1)-minovl, maxshift];
img1c = padarray(img1c, pm, 0, 'post');
img2c = padarray(img2c, pm, 0, 'post');
pm0 = [0, pm(2:end)];
img1c = padarray(img1c, pm0, 0, 'pre');
img2c = padarray(img2c, pm0, 0, 'pre');

img1c = double(img1c); img2c = double(img2c);

cm = getParameter(param, 'mean', true);
if cm
   img1c = img1c - mean(img1c(:));
   img2c = img2c - mean(img2c(:));
end


% weights
pm1 = [1, pm(2:end) + 1];

w1 = zeros(size(img1c));
w1 = imreplace(w1, ones(s1c), pm1);

if isequal(s1c, s2c)
   w2 = w1;
else
   w2 = zeros(size(img2c));
   w2 = imreplace(w2, ones(s2c), pm1);
end

%figure(96); clf; implottiling({mat2gray(w1), mat2gray(img1c); mat2gray(w2) , mat2gray(img2c)});

%fft
w1fft = fftn(w1); 
if isequal(s1,s2)
   w2fft = w1fft;
else
   w2fft = fftn(w2);
end
i1fft = fftn(img1c); i2fft = fftn(img2c);

wssd = ifftn(i1fft .* conj(i2fft));
nrm  = ifftn(w1fft .* conj(w2fft));

%figure(97); clf; implottiling({mat2gray(wssd), mat2gray(nrm), mat2gray(wssd ./ nrm)});

% extract range of interest
wssd = circshift(wssd, pm0);
nrm  = circshift(nrm,  pm0);
nrm(nrm <= 0) = eps;

pm2 = pm; pm2(2:end) = 2 * pm2(2:end);
wssd = imextract(wssd, ones(1, length(pm)), pm2 + 1);
nrm  = imextract(nrm , ones(1, length(pm)), pm2 + 1);

%figure(98); clf; implottiling({mat2gray(wssd), mat2gray(nrm), mat2gray(wssd ./ nrm)});
%size(wssd)

% uable  size is first coordinate and p-m in both directions
cc = wssd ./ nrm;

[quality, pos] = max(cc(:));
pos = imind2sub(size(cc), pos);

% correct for image cutting

shift = zeros(1,dim);
shift(1) = s1(1) - maxovl - 1;
shift = shift + pos - pm0 - 1;

end