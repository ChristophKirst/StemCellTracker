function [shift, quality] = align2ImagesByRMS(img1, img2, varargin)
%
% [shift, quality] = align2ImagesByRMS(img1, img2)
%
% description: 
%     globally align two images using overlap weighted root mean square differences calculated via fft
%
% input:
%     img1,img2    images
%     param        parameter struct with entries
%                  .overlap.min    minimal overlap of the images
%
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units
%     quality      (optional) quality measure of the alignment

param = parseParameter(varargin{:});

% image sizes and dim
s1 = size(img1); s2 = size(img2);
dim = length(s1);
if dim ~= length(s2)
   warning('alignByRMS: image dimension mismatch: %g ~= %g, padding with ones!', length(s1), length(s2));
   % fill smaller dim with ones
   dim = max(dim, length(s2));
   s1 = padright(s1, dim - length(s1), 1);
   s2 = padright(s2, dim - length(s2), 1);
end

% get overlap and regularize
minovl = getParameter(param, 'overlap.min', 1);
if minovl <= 0
   warning('alignByRMS: overlap.min has to be at least 1!');
   minovl = 1;
end
if length(minovl) ~= dim
   pd = minovl(2:end); 
   if isempty(pd)
      pd = minovl(1);
   end
   minovl = padright(minovl, dim, pd);
end
si = min(s1, s2);
for d = 1:dim
   if minovl(d) > si(d)
      minovl(d) = si(d);
   end
end
      

%pad one image if it is not the same size
si = max(s1,s2);
img1 = padarray(img1, si - s1, 0, 'post');
img2 = padarray(img2, si - s2, 0, 'post');

%fprintf('alignByRMS: image sizes: %s, %s, joint: %s\n', var2char(s1), var2char(s2), var2char(si));

% pad border to ensure no overlaps for maximal shifts
pm = (si-minovl);
pmh = ceil(pm/2); 
img1 = padarray(img1, pmh, 0, 'both');
img2 = padarray(img2, pmh, 0, 'both');
img1 = double(img1); img2 = double(img2);

% weights
w1 = zeros(size(img1));
w1 = imreplace(w1, ones(s1), pmh+1);

if isequal(s1, s2)
   w2 = w1;
else
   w2 = zeros(size(img2));
   w2 = imreplace(w2, ones(s2), pmh+1);
end

%figure(96); clf; implottiling({mat2gray(w1), mat2gray(img1); mat2gray(w2) , mat2gray(img2)});

%fft
w1fft = fftn(w1); 
if isequal(s1,s2)
   w2fft = w1fft;
else
   w2fft = fftn(w2);
end
i1fft = fftn(img1); i2fft = fftn(img2);
s1fft = fftn(img1 .* img1);  s2fft = fftn(img2 .* img2);
wssd = w1fft .* conj(s2fft) + s1fft .* conj(w2fft) - 2 * i1fft .* conj(i2fft);
%wssd = i1fft .* conj(i2fft); 

wssd = ifftn(wssd);
nrm  = ifftn(w1fft .* conj(w2fft));

%figure(97); clf; implottiling({mat2gray(wssd), mat2gray(nrm), mat2gray(wssd ./ nrm)});

% extract overlapping range 
pm1 = s1 - minovl; pm2 = s2 - minovl;
wssd = circshift(wssd, pm2);
nrm  = circshift(nrm,  pm2);
nrm(nrm <= 0) = eps;

%figure(98); clf; implottiling({mat2gray(wssd), mat2gray(nrm), mat2gray(wssd ./ nrm)});
%size(wssd)

wssd = imextract(wssd, ones(1, length(pm)), pm1 + pm2 + 1);
nrm  = imextract(nrm , ones(1, length(pm)), pm1 + pm2 + 1);

%figure(99); clf; implottiling({mat2gray(wssd), mat2gray(nrm), mat2gray(wssd ./ nrm)});
%size(wssd)

% uable  size is first coordinate and p-m in both directions
cc = wssd ./ nrm;

[quality, pos] = min(cc(:));
quality = - quality;
pos = imind2sub(size(cc), pos);

% center 
%shift = {zeros(1, length(pos)), pos - pm2 - 1};
shift = pos - pm2 - 1;

end