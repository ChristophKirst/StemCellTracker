function [imglab, varargout] = fucci_label(img, varargin)

if nargin > 1
   verbose = varargin{1};
else
   verbose = false;
end

%% mask

%imgmask = img > 0.125;
%imgmask = img > 0.25;
imgflat = img - imopen(img, strel('disk', 50));
imgmask = imgflat > 0.05;

imgmask = imopen(imgmask, strel('disk', 3)) > 0;
%imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 50) > 0;

if verbose
   %max(img(:))
   
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implot(imoverlaylabel(imgflat, double(imgmask)))
end
   

%%
imgf =imgflat;
imgf = filterLoG(max(imgf(:)) - imgf, [38,38], []);
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
imgmax = imdilate(imgmax, strel('disk', 6));


% constrain to maxima within mask
imgmax = immask(imgmax, imgmask);

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

% filter small peaks

imglab = postProcessSegments(imglab, 'volume.min', 100);

if verbose
   figure(23); clf;
   set(gcf, 'Name', ['Labeling'])
   implot(imoverlaylabel(img, imglab))
end

%%
if nargout > 1
   varargout{1} = imgmask;
end