function [imgpost, stats] = segment2dDilute(img, verbose, figure_offset)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% segmentation of diluted cells for tracking %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
   verbose = false;
end
if nargin < 3
   figure_offset = 0;
end

img = mat2gray(img);


if verbose
   figure(1 + figure_offset); clf
   set(gcf, 'Name', 'Image');
   colormap(gray)
   implot(img)
   
   figure(2 + figure_offset); clf
   set(gcf, 'Name', 'Histogram');
   subplot(1,2,1);
   hist(img(:), 256)
   subplot(1,2,2);
   hist(log(img(:))+eps, 256)
end

%%

th = 2^-2.8;
imgth = img;
imgth(imgth < th) = th;
imgth = mat2gray(imgth);

imgmask = img > th;
%imgmask = imclose(imgmask, strel('disk', 1));
imgmask = imopen(imgmask, strel('disk', 5));

imglab = bwlabeln(imgmask);
imglab = postProcessSegments(imglab, setParameter('volume.min', 150, 'relabel' , false));
imgmask = imglab > 0;

if verbose
   figure(3 + figure_offset); 
   set(gcf, 'Name', 'Thresholded Image');
   colormap gray
   implot(imgth)
   
   
   figure(4 + figure_offset)
   set(gcf, 'Name', 'Mask');
   colormap gray
   implottiling({img, imoverlay(imgth, imgmask, 'r',  false)})
end

%% Seeding
imgs = imgth;
imgs = medianFilter(imgs);

imgs = diskFilter(imgs, [15, 15], 1, 1, 0);
imgs(imgs< 0) = 0;
imgs = mat2gray(imgs);

imgmax = imextendedmax(imgs, 0.01);


if verbose 
   figure(5 + figure_offset); clf
   set(gcf, 'Name', 'Seeds');
   colormap gray
   implottiling({imoverlay(imgth, imgmax, 'r', true), imoverlay(imgs, imgmax, 'r', true)})
end

%% Watershed

imgw = imgth;

imgmin = imimposemin(max(imgw(:)) - imgw, imgmax);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);
imgseg = imrelabel(imgseg);

if verbose
   figure(6 + figure_offset); clf
   set(gcf, 'Name', 'Watershed');
   colormap jet
   implottiling({imoverlay(imgw, imgmax), imcolorize(imgseg), imoverlaylabel(img, imgseg, true)});
end

%% Postprocess and create intial statistics

param = setParameter('volume.min',    0,...     % minimal volume to keep (0)
                     'volume.max',    inf,...    % maximal volume to keep (Inf)
                     'intensity.min', -inf, ...  % minimal mean intensity to keep (-Inf)
                     'intensity.max', inf, ...   % maximal mean intensity to keep (Inf)
                     'boundaries',    false, ... % clear objects on x,y boundaries (false)
                     'fillholes',     true,...   % fill holes in each z slice after processing segments (true)
                     'relabel',       true);    % relabel from 1:nlabelnew (true)

[imgpost, stats] = postProcessSegments(imgseg, param);


if verbose 
   figure(7 + figure_offset)
   set(gcf, 'Name', 'Final');
   colormap jet
   implottiling({imoverlay(img, imgmax), imcolorize(imgpost), imoverlaylabel(img, imgpost, true)});
   
   fprintf('seeds: %g\n', max(imgseg(:)));
   fprintf('after postprocessing: %g\n', max(imgpost(:)));
end

end