%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [imglab, varargout] = fucci_label(img, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
end

%% mask

%imgmask = img > 0.125;
%imgmask = img > 0.25;
imgmask = img > 0.1;
imgmask = imopen(imgmask, strel('disk', 3));
imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 500) > 0;

if verbose
   %max(img(:))
   
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implot(imoverlaylabel(img, double(imgmask)))
end
   
imgf =immask(img, imgmask);

%% edge detection

imgf2 = medianFilter(imgf, [5,5]);

imge = edge(imgf2, 'sobel')

figure(3); clf; imcolormap('gray');
implottiling({imgf, imge})





%% peak detection
imgf = logFilter(max(imgf(:)) - imgf, [150,150], []);
%imgf = mat2gray(imgf);

%{min(imgf(:)), max(imgf(:))}

% h-max detection (only local maxima with height > hmax are considered as maxima
%imgmax = imextendedmax(imgf,  0.00005);
imgmax = imextendedmax(imgf,  0.000005);

%local max
%imgmax = imregionalmax(imgf);

% constrain to maxima within mask
imgmax = immask(imgmax, imgmask);

% Combine nearby points
%imgmax = imdilate(imgmax, strel('disk', 4));
imgmax = imdilate(imgmax, strel('disk', 5));


% fill holes - combination of nearby points can lead to holes
imgmax = imfill(imgmax,'holes');
% shrink to single points - extended maxima usually give better segmentation results
% imgmax = bwmorph(imgmax,'shrink',inf);           

% plot the results.
if verbose  
   figure(22); clf
   set(gcf, 'Name', ['Seeding'])
   implottiling({imoverlay(img, imgmax), imoverlay(imgf, imgmax)});
end


%% Label

imglab = bwlabeln(imgmax);

if verbose
   figure(23); clf;
   set(gcf, 'Name', ['Labeling'])
   implot(imoverlaylabel(img, imglab))
end

%%
if nargout > 1
   varargout{1} = imgmask;
end