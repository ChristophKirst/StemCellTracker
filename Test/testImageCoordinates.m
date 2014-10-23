%% Note on Image Coordinates
%
% Synposis: spatial coordinates (x,y,z)
%           image coordinates (p,q,l) [= (x,y,z) in Imtools if not specified otherwise]
%           image sizes (w,h,d)
%           intensity at (x,y,z) is img(x,y,z) 
%
% Matlab represents arrays starting with columns.
% Normal coordinate specification x,y,z is assumed in the ImTools
% i.e. for an image img, the pixel img(p,q) is at (x,y) = (p,q)
% and img(p,q,l) is at (x,y,z) = (p,q,l)
%
% ImTools uses pixel coordinates (p,q,l) consistently in all its functions
% In particular implot will be consistetn with adding lines in p,q,l coordinates
% A 2d matrix img will be plotted effectively as img', ad 3d matrix as permute(img, [2,1,3])
% to make this work.
%
% Pixel Coordinates are denoted as: p,q,l, spatial coordinates as x,y,z 
% (which coincide with pixel coordinates in the ImTools
% but not necessarily with other matlab functions such as regionprops.)
% All internal ImTools routines correct for the matlab crazyness if necessary
% and give consistent results.
%
% Sizes in different image dimensinos are denoted as w, h, d (for width, height, depth)
%

%% Example: Consistent coordinates using ImTools

test = zeros(40,30);
test = imreplace(test, fspecial2('disk', [20,10],0), [5, 10]);
test = mat2gray(test);

p = 33; q = 21;
test(p, q) = 1;
implot(test)
line([5, p], [10, q])


%% Example: Matlab is unintuitive

test = zeros(40,30);
test = imreplace(test, fspecial2('disk', [20,10],0), [5, 10]);
test = mat2gray(test);

p = 33; q = 21;
test(p, q) = 1;
imshow(test)
line([5, p], [10, q])


%% Examples: improfile uses matlab x,y coordinates that are inconsistent with pixel coordinates

test = 1:30;
test =repmat(test', 1, 20);
test = repmat((1:20),30, 1) .* test;
[sp, sq] = size(test)

x = [6 17]; % spatial coordinates
y = [2 23]; % spatial coordinates

pq = [y; x]'  % corresponding pixel coordinates
p = y;
q = x;

mask =imcoords2mask([sp,sq], pq);

figure(40)
subplot(1,2,1)

imshow(imoverlay(mat2gray(test), mask))
line(x,y)   % uses space coordinates
title(sprintf('image size is: hxw = %gx%g', sp, sq))

pp = improfile(test, x, y);  % uses space coordinates 
ep = test(sub2ind(size(test), p, q)); % pixel indices are in h,w coordinates !

subplot(1,2,2)
plot(pp)
title(sprintf('start end intensities: %g %g', ep(1), ep(2)))

figure(41)
clf
implot(imoverlay(mat2gray(test), mask))
set(gca,'YDir','Reverse')
line(p,q)

