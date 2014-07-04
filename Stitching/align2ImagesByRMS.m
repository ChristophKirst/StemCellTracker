function shift = align2ImagesByRMS(img1, img2, varargin)
%
% shift = align2ImagesByRMS(img1, img2)
%
% description: 
%     align 2 images using cross correlation calculated via fft
%
% input:
%     img1,img2    images
%     param        parameter struct with entries
%                  .overlap.min   minimal overlap of images (50)
%                  .overlap,max   maximal overlap considered (150)
%                  .shift.max     maximal shift in secondary directions ([10, 10])
%                  .shift.min     minimal shift in secondary directions ([-10, -10])
%
% output:
%     shift        the shift between origin of img1 to origin of img2 in pixel coordinates and pixel units

param = parseParameter(varargin{:});

minov = getParameter(param, 'overlap.min', 50);

%pad one image if it is not the same size
s1 = size(img1); s2 = size(img2);
si = max(s1,s2);

img1 = padarray(img1, si - s1, 0, 'post');
img2 = padarray(img2, si - s2, 0, 'post');

img1  = padarray(img1, si, 0, 'both');
img2  = padarray(img2, si, 0, 'both');
img1 = double(img1); img2 = double(img2);


% weights
w = zeros(size(img1));
w = imreplace(w, ones(si), si+1);

%fft
wfft = fftn(w);
i1fft = fftn(img1); i2fft = fftn(img2);
s1fft = fftn(img1 .* img1);  s2fft = fftn(img2 .* img2);
wssd = wfft .* conj(s2fft) + s1fft .* conj(wfft) - 2 * i1fft .* conj(i2fft);
%wssd = i1fft .* conj(i2fft);
%max(abs(wssd(:)))
%min(abs(wssd(:)))

wssd = ifftshift(ifftn(wssd));
nrm = ifftshift(ifftn(wfft .* conj(wfft)));

sd2 = fix(si/2);
cc = imextract(wssd, sd2 + 1 + minov, sd2 + 2 * si - 1 - minov) ./ imextract(nrm, sd2 + 1 + minov, sd2 + 2 * si - 1 - minov);

[~, pos] = min(cc(:));
pos = imind2sub(size(cc), pos);

shift =  pos - si + minov - 1;

end