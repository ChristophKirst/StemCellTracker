%%%%%%%%%%%%%%%%%%%%%%%
%%% Find Thresholds %%%
%%%%%%%%%%%%%%%%%%%%%%%

% this script is intened to help find a good thresdhold for masking


%% direct thresholding methods (with optional prefiltering)

%param.filter.ksize = 3;
%imgf = medianFilter(img, param.filter.ksize);
imgf = img;

thOtsu = thresholdOtsu(imgf)
thEntropy = thresholdEntropy(imgf)
thMutualEntropy = thresholdMutualEntropy(imgf)
thMOG = thresholdMixtureOfGaussians(imgf, 0.5)



%% determine threshold using the histogram on logarithmic intensities
imgvalslog = log2(img(:)+eps);
imgvalslog(imgvalslog < -15) = -15; % bound lower values
imgvalslog(imgvalslog > 0) = 0;     % bound upper values

thLogMOG = 2^thresholdMixtureOfGaussians(imgvalslog, 0.5)  % this usually takes long


%% thresholding using local maxima statistics

imgvalsmax = img(imregionalmax(img));
%imgvalsmax = img(imextendedmax(img, 0.01));
imgvalslogmax = log2(imgvalsmax);
 
thMaxMOG = 2^thresholdMixtureOfGaussians(imgvalslogmax, 0.5)

%% select a threshold and create mask and thresholded image
th = thMaxMOG;
th = 0.15;

imgth = img;
imgth(imgth < th) = 0;
imgmask = imgth > 0;

% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 2));    % larger disk size removes larger fragments% close
%imgmask = imclose(imgmask, strel('disk', 2));  % larger disksize closes more
% dilate
%imgmask = imdilate(imgmask, strel('disk', 2)); % increase mask region 
% erode
%imgmask = imerode(imgmask, strel('disk', 2));  % decreasde mask region
% fill holes
%imgmask = imfill(imgmask,'holes');

% plot results

figure_offset = 100;

prt = '\n\nthresholds:\n===========\n';
thnames = {'thOtsu', 'thEntropy', 'thMutualEntropy', 'thMOG', 'thLogOtsu', 'thLogEntropy', 'thLogMutualEntropy', 'thLogMOG',  'thLogOtsu', ...
           'thMaxOtsu', 'thMaxEntropy', 'thMaxMutualEntropy', 'thMaxMOG'};

for t = 1:length(thnames)
   if exist(thnames{t}, 'var')
      prt = [prt '\n' sprintf(['%18s = %7.5f'], thnames{t}, eval(thnames{t}))];
   end
end
prt = [prt '\n----------------------------\n'];
prt = [prt sprintf('%18s = %7.5f', 'th', th)];
fprintf([prt '\n']);


figure(1 + figure_offset)
set(gcf, 'Name', ['Thresholding: ' filename])
implottiling(imgth)

figure(2 + figure_offset)
set(gcf, 'Name', ['Mask: ' filename])
implottiling(imoverlay(img, imgmask));

figure(3 + figure_offset)
set(gcf, 'Name', ['Thresholding Histograms: ' filename])

subplot(2,4,1);
hist(img(:), 256);
title('raw intensities')
subplot(2,4,5);
plot(sort(img(:))) % x-axis is effectively number of pixels <= ordinate in following plots

if exist('imgvalslog', 'var')
   subplot(2,4,2);
   hist(imgvalslog, 256);
   title('log intensities');
   subplot(2,4,6);
   plot(sort(imgvalslog(:)))
end
if exist('imgvalsmax', 'var')
   subplot(2,4,3);
   hist(imgvalsmax, 256);
   title('local max intensities');
   subplot(2,4,7);
   plot(sort(imgvalsmax(:)))
end
if exist('imgvalslogmax', 'var')
   subplot(2,4,4);
   hist(imgvalslogmax, 256);
   title('local log max intensities');
   subplot(2,4,8);
   plot(sort(imgvalslogmax(:)))
end

