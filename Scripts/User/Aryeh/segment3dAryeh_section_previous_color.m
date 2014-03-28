%function imgseg = segment3d(img, verbose)
%
% 3d segmentation for geometrically confined stemm cell colonies
%

%% Initialize Segmenter
if false
    %%
    close all
    clear all
    clc
    set(0, 'DefaultFigurePosition', [1   705   560   420]);
end
initialize

%% Init

if nargin < 2
    verbose = false;
end

useimaris = false;
savefile = '';
verbose = false;

initialize();

%%
if false
    %% Init for Testing
    bfinitialize();

    %%
    if useimaris
       %%
       imarisinitialize();
    end

    %% Data form disk

    %filename = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2.lif';
    filename = '\\tracking-pc.rockefeller.edu\DATA\Aryeh\hESC\Cytoo_IF\Expt24_Snail_Bra_Sox2_timecourse\140305_RUES2_36hBMP4_Bra_Snail_Sox2.lif';

    %lifdata = imread_bf(filename, struct('series', 21, 'channel', 1));
    %lifdata = imread_bf(filename, struct('series', 21, 'time', 1, 'channel', 1, 'y', 1000 + [0, 512]));
    %img = lifdata(:,:,:,1,1);
    %clear lifdata

    %% samle tile
    xr = [1, 512];  % use [] for all
    yr = [1, 512];
    cr = 1;
    se = 2;
    ti = 1;
    img = imread_bf(filename, struct('series', se, 'time', ti, 'channel', cr, 'x', xr, 'y', yr));
    img = imzreverse(squeeze(img));


    %%  figure aryeh
    xr = [275, 2100];  % use [] for all
    yr = [825, 1175];
    cr = 1;
    se = 21;
    ti = 1;
    img = imread_bf(filename, struct('series', se, 'time', ti, 'channel', cr, 'x', xr, 'y', yr));
    img = imzreverse(squeeze(img));


    %% Data from imaris

    if useimaris
       %%
       imarisstart();

       %%
       img = imarisget('Volume', 0,0);

       %%

       imarissetvolume(uint8(img))

    end


else

    filename = '';
end


%% Plot low resolution
if verbose
    %%
    figure(1)
    clf
    downsamplexy = 5;
    imgres = img(1:downsamplexy:end, 1:downsamplexy:end, 1:1:end);
    implot3d(mat2gray(imgres))
end



%% initialize / prefiltering

imgd = double(img);
imgd = mat2gray(imgd);

%imglogvals = log2(imgd(:)+eps);
%imglogvals(imglogvals < -5) = -5;
%imglogvals(imglogvals > 0) = 0;

if verbose
    figure(2)
    clf
    set(gcf, 'Name', ['Raw Stack: ' filename ' channel: 1']);
    implot3d(imgd);

    %figure(2)
    %subplot(1,2,1);
    %hist(imglogvals, 256)
    %subplot(1,2,2);
    %hist(imgd(:), 256);
end


%% filter
param.filter.median.ksize = [3, 3, 3]; %3;
imgmed = mat2gray(medianFilter(imgd, param.filter.median.ksize));
%imgmed = meanShiftFilter(imgmed, 4, 0.1);


imgmedf = imgmed;
imgmedf(imgmed > 0.6) = 0.6;
imgmedf = mat2gray(imgmedf);

if verbose
    %%
    figure(3)
    clf
    set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
    implot3d(imgmed);

    figure(4)
    clf
    set(gcf, 'Name', ['Filtered Stack: ' filename ' channel: 1']);
    implottiling(imgmedf);


    %figure(5)
    %hist(imgmed(:), 256)
end





%% thresholding / masking

%imglogvals = log2(imgmedf + eps);
%imglogvals(imglogvals < -6) = -6;

%th = 2^thresholdMixtureOfGaussians(imglogvals, 0.9)
%th = thresholdMixtureOfGaussians(imgmedf, 0.5)
%th = thresholdEntropy(imgmed)
%th = thresholdMutualEntropy(imgmed)
%th = thresholdOtsu(imgmed);
th = 0.1;

imgth = imgmedf;
imgth(imgth < th) = 0;

imgmask = imgth > 0;
% remove small fragments and small isolated minima with imopen
imgmask = imopen(imgmask, strel('disk', 3));


if verbose
    %%
    figure(4)
    clf
    set(gcf, 'Name', ['Thresholded Stack: ' filename ' channel: 1']);
    implot3d(imgth);

    %figure(5)
    %clf
    %set(gcf, 'Name', ['Mask: ' filename ' channel: 1']);
    %imsurfaceplot3d(imgmask);   lnew = lnew + 1;

    %ijplot3d(imgth, 'PixelDepth', 5)



    %%
    figure(5)

    subplot(1,2,1)
    hist(imglogvals(:), 256/8)

    subplot(1,2,2);
    hist(imgmedf(:), 256/8);
end


%% Filtering / Seeding

% gaussian smoothing
%imgdg = gaussianFilter(imgd,3,10);

% median filter  ?? also computed above
% imgf = medianFilter(imgdg, 3);

% mean shift
%imgf = meanShiftFilter(imgd, 3, 0.1);

% disk [ksize or outer-box, width of ring, wt-inner disk, wt-ring]
%imgf = diskFilter(imgmed, [20, 20, 4], 1, 1, -1);

% Laplacian of Gaussians

zsize = size(imgmedf,3);
imgf = mat2gray(imgmedf);
%for z = zsize:-1:1
%   imgf(:,:,z) = imgf(:,:,z) - 0.0 * mat2gray(imgradient(imgf(:,:,z)));
%end
imgf = mat2gray(imgf);


if verbose
    figure(17)
    clf
    imsubplot(1,2,1);
    implot3d(imgf);
end

%param.filter.logsize = [10,10,10];
%imgf = logFilter(max(imgf(:)) - imgf, param.filter.logsize, [], 0);

imgf = dogFilter(imgf, [15, 15, 7], [] ,[], 0);


if verbose
    imsubplot(1,2,2);
    imgs = imgf;
    imgs(imgs< 0) = 0;
    implot3d(mat2gray(imgs));
end


%imgf = sphereFilter(imgf, [7,7,7], [], 0);

%imgf = diskFilter(imgf, [11,11,9], 1, 1, -1, 0);
%imgf(imgf < 0 ) = 0;
imgf = imgf - min(imgf(:));
imgf = mat2gray(imgf);




%%%
% find maxima using h-max detection
param.filter.hmax = 0.01;  %0.02;
%imgmax = imregionalmax(imgf);
imgmax = imextendedmax(mat2gray(imgf), param.filter.hmax);
imgmax = imdilate(imgmax, strel('disk',2));

if verbose
    %%
    figure(21)
    clf
    implottiling(imoverlaylabel(mat2gray(imgmedf), bwlabeln(imgmax)));

    figure(22)
    clf
    implottiling(imoverlaylabel(mat2gray(imgf), bwlabeln(imgmax)));

end


%% Watershed Segmentation on Image

%dilating the maxima can help
%imgmaxd = imdilate(imgmax, strel('disk', 1));
imgmaxd = imfill(imgmax, 'holes');
%imgmaxd = imgmax;

% watershed
imgmin = imimposemin(max(imgmedf(:)) - imgmedf, imgmaxd);
imgws = watershed(imgmin);
imgseg = immask(imgws, imgmask);

%propagation
%imgseg = segmentByPropagation(imgmedf, bwlabeln(imgmaxd), imgmask);



% %% Clean up segmentation and alternative diagnositcs
%
imgseg = imlabelseparate(imgseg);
stats = regionprops(imgseg, 'Area', 'PixelIdxList');

min_vol = 100;
keep = [stats.Area] >= min_vol;
for i = find(~keep)
    imgseg(stats(i).PixelIdxList) = 0;
end
imgseg = imfill(imgseg, 'holes');
imgseg = imrelabel(imgseg);

%
% [bndry, lbl] = bwboundaries(lbl > 0, 8, 'noholes');
% stats = regionprops(logical(lbl), 'Area', 'Centroid', 'PixelIdxList');
% imgmask = lbl > 0;

if verbose
    %%
    figure(11)
    clf
    %implot3d(imgseg);
    imsurfaceplot3d(imgseg);
end


%%
if verbose
    %%
    figure(23)
    clf
    implotlabeloutline(img, imgseg);
end


%%
if true
    ijplot3d(imcolorize(imgseg) + gray2rgb(mat2gray(img)), 'PixelDepth', 5);
end

%%
size(imlabel(imgseg))
max(imgseg(:))








%% Calculate Surfaces

[surf, fac, norm] = imsurface(imgseg, 'all');
surfaces = {surf, fac, norm};


%% Push Surfaces to imaris

if useimaris
    %%
    imarisstart
    
    
   %%
   %nset = 3;
   nset = size(surf);
   sfset = surf(1:nset);
   fcset = fac(1:nset);
   nmset = norm(1:nset);

   imarissetsurface('Segmentation', sfset, fcset, nmset, 0, [1,1,1]);
end




%% Calcualte Statistics 

stats = cell(4,1);
stats{1} = statisticsSegments(double(img), imgseg);

%%
imgc{1} = imgcl;

%clip = {[42 165], [20, 93], [34, 71], [42, 65]};
clip = {[0, 255], [24 70], [14, 58], [21, 71]};

for ch = 2:4

   imgc{ch} = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
   imgc{ch} = imzreverse(squeeze(imgc{ch}));
   imgc{ch} = imisostack(imgc{ch}, 1, 'min');
   imgc{ch} = imclip(imgc{ch}, clip{ch});
   imgc{ch} = mat2gray(double(imgc{ch}));
   
   stats{ch} = statisticsSegments(imgc{ch}, imgseg);
end

cdapi = 1;
cbra = 2;
csnail = 3;
csox2 = 4;

%% save this stuff

if ~isempty(savefile)
   %%
   %save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_imaris.mat', 'imgseg')

   %save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris_figure.mat', 'surfaces')

   save('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_imgseg_imaris_section.mat', 'imgseg')
   save('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris_section.mat', 'surfaces')
   
   %%
   save(savefile, 'img', 'imgseg', 'stats', 'surfaces')
   
   
   %%
   save('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_imgseg_imaris_matlab_session_previous.mat')
end

%% load

load('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_imgseg_imaris_matlab_session_previous.mat')





%%
imgc{1} = img;

%clip = {[42 165], [20, 93], [34, 71], [42, 65]};
clip = {[0, 255], [24 70], [14, 58], [21, 71]};

for ch = 2:4

   imgc{ch} = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
   imgc{ch} = imzreverse(squeeze(imgc{ch}));
   imgc{ch} = imisostack(imgc{ch}, 1, 'min');
   imgc{ch} = imclip(imgc{ch}, clip{ch});
   imgc{ch} = mat2gray(double(imgc{ch}));
   
   
   %set the bottom three layers of sox2 to zero
   if ch==4
      img2 = imgc{ch};
      img2(:,:,1:4) = 0;
      imgc{ch} = img2;
   end
   
   stats{ch} = statisticsSegments(imgc{ch}, imgseg);
end

cdapi = 1;
cbra = 2;
csnail = 3;
csox2 = 4;


%% select subset of surfaces of interset





%% Visualization / Plotting Statistics etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nset = size(surf);
sfset = surf(1:nset);
fcset = fac(1:nset);
nmset = norm(1:nset);

imarissetsurface('Segmentation', sfset, fcset, nmset, 0, [1,1,1]);


%% Color code

% 2 = bra = green
% 3 = snail = red
% 4 = sox2 = blue

val = 'MeanIntensity';
b = [stats{2}.(val)]';
r = [stats{3}.(val)]';
g = [stats{4}.(val)]';

% mm = 40; mx = 126;
% r(r < mm) = mm; r(r> mx) = mx;
% 
% mm = 50; mx = 107;
% g(g < mm) = mm; g(g> mx) = mx;
% 
% mm = 20; mx = 92;
% b(b < mm) = mm; b(b> mx) = mx;

g = mat2gray(g); r = mat2gray(r); b = mat2gray(b);
%g = g/256; r = r/256; b = b/256;
%g = g/max(g(:)); b = b/max(b(:)); r = r/max(r(:)); 
%r = r/256;
%
%b = zeros(size(r,1),1);
%r = zeros(size(r,1),1);
%g = ones(size(r,1),1);

%nsurfaces = length(stats{1});
%r = (0:nsurfaces-1)' / (nsurfaces-1);
%g = 0 * r; b = 0 * r;

%make it more light
%fbirght = 2;
%r = fbirght*r; g = fbirght*g; b = fbirght*b;

%fshift = 0.5; %-> pastel effect
%r = r + fshift; g = g + fshift; b = b + + fshift;
%r(r>1) = 1; g(g>1) = 1; b(b>1) = 1;

%rgb = [r,g,b];
%grayid = (sqrt(sum(rgb .* rgb, 2)) < 0.01);
%r(grayid) = 0.75; g(grayid) = 0.75; b(grayid) = 0.75;


if true
%%
figure(13)
subplot(1,3,1)
hist(r,256)
subplot(1,3,2)
hist(g,256)
subplot(1,3,3)
hist(b,256)

end


% label accrording to fate

%snail > th -> red
rid = r>0.05;
r(~rid) = 0;
r(rid) = 1; 


bid = b > 0.015;
b(~bid) = 0;
b(bid) = 1;

gid = g > 0.15;
g(~gid) = 0;
g(gid) = 1;
g(rid) = 0;
b(gid) = 0;

zid = and(and(r == 0, g == 0), b == 0);
g(zid) = 1;r(zid) = 1;b(zid) = 1;


c = imarisrgb2statistics(r,g,b);
imarissetstatistics('rgbcolor20', c);
 




%% Plotting Statistics as Colored Surfaces

if verbose
   %%
   figure(11)
   clf
   %implot3d(imgseg);
   %cdata = [stats.Centroid];
   %cdata = cdata(3,:);
   
   %cdata = [stats.MeanIntensity];
   cdata = [stats{3}.MeanIntensity];

   %cdata = [stats{2}.Volume];
   
  
   param.color.data = cdata;
   imsurfaceplot3d(imgseg, param)
end

%% Plotting Statistics as Scatter Plots
if verbose
   
   xmeasure = 'Volume';
   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   cccol = {'b', 'g', 'r', 'k'};
   
   for cc = 1:4
      figure(60+cc)
      clf
      set(gcf, 'Name', ['channel: ' num2str(cc)]);
      for i = 1:4
         ax(i) = subplot(2,2,i);
         plot([stats{cc}.(xmeasure)], [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
         xlabel(xmeasure); ylabel(ymeasures{i});
      end
      linkaxes(ax, 'x')
   end  
end


%%
if verbose
   
   xmeasure = 'z';
   xdata = [stats{1}.Centroid];
   xdata = xdata(3,:);

   ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
   cccol = {'b', 'g', 'r', 'k'};

   
   for cc = 1:4
      figure(70+cc)
      clf
      set(gcf, 'Name', ['channel: ' num2str(cc)]);
      for i = 1:4
         ax(i) = subplot(2,2,i);
         plot(xdata, [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
         xlabel(xmeasure); ylabel(ymeasures{i});
      end
      linkaxes(ax, 'x')
   end  
   
end


%% Plot Statistics in Imaris

%zsl = 3:size(img,3);
zsl = 1:size(img,3);

sz= size(img(:,:,zsl));
nch = 4;

dset = imarissetdataset('uint8', sz(1), sz(2), sz(3), nch, 1);

imarissetvolume(dset, uint8(img(:,:,zsl)),0)

for ch = 2:4

   img2 = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
   img2 = imzreverse(squeeze(img2));
   
   %img2 = imgc{ch};
   if ch==4
      img2(:,:,1:3) = 0;
   end
   
   imarissetvolume(dset, uint8(img2(:,:,zsl)),ch-1)
end


%%
vRGBA = imarisrgb2color(1, 1, 1);
dset.SetChannelColorRGBA(0, vRGBA);

vRGBA = imarisrgb2color(0, 0, 1);
dset.SetChannelColorRGBA(1, vRGBA);

vRGBA = imarisrgb2color(1,0,0);
dset.SetChannelColorRGBA(2, vRGBA);

vRGBA = imarisrgb2color(0, 1, 0);
dset.SetChannelColorRGBA(3, vRGBA);



%% set colors

vRGBA = imarisrgb2color(0, 0, 1);
dset.SetChannelColorRGBA(0, vRGBA);

vRGBA = imarisrgb2color(1, 1, 1);
dset.SetChannelColorRGBA(1, vRGBA);

vRGBA = imarisrgb2color(1,0,0);
dset.SetChannelColorRGBA(2, vRGBA);

vRGBA = imarisrgb2color(0, 1, 0);
dset.SetChannelColorRGBA(3, vRGBA);

%% threshold in sox2 / snail intensity

figure(120)
subplot(1,2,1)
hist([stats{3}.MedianIntensity], 255)
xlabel('Snail')

subplot(1,2,2)
hist([stats{4}.MedianIntensity], 255)
xlabel('Sox2')


%%

figure(120)
subplot(1,2,1)
hist([stats{3}.MeanIntensity], 255)
xlabel('Snail')

subplot(1,2,2)
hist([stats{4}.MeanIntensity], 255)
xlabel('Sox2')



%%
%nset = 3;
nset = size(surf);
sfset = surf(1:nset);
fcset = fac(1:nset);
nmset = norm(1:nset);

imarissetsurface('Segmentation', sfset, fcset, nmset, [1, 1, 0.5]);


%% 

thsnail = 8.61;
thsox2 = 23.3;


snail = [stats{3}.MeanIntensity];
sox2   = [stats{4}.MeanIntensity];
ncells = length(snail);
snailsox2 = zeros(1, ncells);
snailsox2(sox2 > thsox2) = 2;
snailsox2(snail > thsnail) = 3;
snailsox2(and(snail > thsnail, sox2 > thsox2)) = 1;

imarissetstatistics('snailsox', snailsox2);


%%

% create separate surfaces for different 

ids = find(and(sox2 >= thsox2, snail < thsnail));

nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('Sox2', sfset, fcset, nmset);

%%

ids = find(and(snail >= thsnail, sox2 < thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('Snail', sfset, fcset, nmset);


%%
ids = find(and(snail >= thsnail, sox2 >= thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('SnailSox2', sfset, fcset, nmset);

%%
ids = find(and(snail < thsnail, sox2 < thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('None', sfset, fcset, nmset);



%end



%%

save('Z:\aryeh_sox2_snail.mat')



%%



% set segmented colored volume:

nsurfaces = length(stats);

rvol = imgseg;
gvol = imgseg;
bvol = imgseg;

for i = 1:nsurfaces
   rvol(stats{i}.PixelIdxList) = r(i);
   gvol(stats{i}.PixelIdxList) = g(i);
   bvol(stats{i}.PixelIdxList) = g(i);
end

rvol = imisostack(rvol,-1, 'mean');
gvol = imisostack(gvol,-1, 'mean');
bvol = imisostack(bvol,-1, 'mean');


% Create 3 new channels

idatset = imarisgetdataset();
nchannels = idatset.GetSizeC;
vDataSet.SetSizeC(4+3);            
vDataSet.SetChannelName(4, 'Segmentation Red');
vDataSet.SetChannelName(5, 'Segmentation Green');
vDataSet.SetChannelName(5, 'Segmentation Green');

































%% Calcualte Statistics

stats = cell(4,1);
stats{1} = statisticsSegments(double(img), imgseg);

%%
for ch = 2:4

    img2 = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
    img2 = imzreverse(squeeze(img2));

    stats{ch} = statisticsSegments(img2, imgseg);
end








%% Calculate Surfaces

[surf, fac, norm] = imsurface(imgseg, 'all');
surfaces = {surf, fac, norm};

%% save this stuff

if ~isempty(savefile)
    %%
%save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_segmetation_imaris.mat', 'imgseg')

%save('./Test/Images/Develop/Aryeh/140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris_figure.mat', 'surfaces')

save('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_imgseg_imaris_section.mat', 'imgseg')
save('Z:\140305_RUES2_36hBMP4_Bra_Snail_Sox2_surfaces_imaris_section.mat', 'surfaces')

    %%
    save(savefile, 'img', 'imgseg', 'stats', 'surfaces')
    
    
    
    
    
end





























%% Visualization / Plotting Statistics etc


%% Push Surfaces to imaris

if useimaris
    %%
    %nset = 3;
    nset = size(surf);
    sfset = surf(1:nset);
    fcset = fac(1:nset);
    nmset = norm(1:nset);

    imarissetsurface('Segmentation', sfset, fcset, nmset);
end


%% Plotting Statistics as Colored Surfaces

if verbose
    %%
    figure(11)
    clf
    %implot3d(imgseg);
    %cdata = [stats.Centroid];
    %cdata = cdata(3,:);

    %cdata = [stats.MeanIntensity];
    cdata = [stats{3}.MeanIntensity];

    %cdata = [stats{2}.Volume];


    param.color.data = cdata;
    imsurfaceplot3d(imgseg, param)
end

%% Plotting Statistics as Scatter Plots
if verbose

    xmeasure = 'Volume';
    ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 
'MedianIntensity'};
    cccol = {'b', 'g', 'r', 'k'};

    for cc = 1:4
       figure(60+cc)
       clf
       set(gcf, 'Name', ['channel: ' num2str(cc)]);
       for i = 1:4
          ax(i) = subplot(2,2,i);
          plot([stats{cc}.(xmeasure)], [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
          xlabel(xmeasure); ylabel(ymeasures{i});
       end
       linkaxes(ax, 'x')
    end
end


%%
if verbose

    xmeasure = 'z';
    xdata = [stats{1}.Centroid];
    xdata = xdata(3,:);

    ymeasures = {'MeanIntensity', 'MinIntensity','MaxIntensity', 'MedianIntensity'};
    cccol = {'b', 'g', 'r', 'k'};


    for cc = 1:4
       figure(70+cc)
       clf
       set(gcf, 'Name', ['channel: ' num2str(cc)]);
       for i = 1:4
          ax(i) = subplot(2,2,i);
          plot(xdata, [stats{cc}.(ymeasures{i})], ['*' cccol{cc}]);
          xlabel(xmeasure); ylabel(ymeasures{i});
       end
       linkaxes(ax, 'x')
    end

end


%% Plot Statistics in Imaris

zsl = 3:size(img,3);

sz= size(img(:,:,zsl));
nch = 4;

dset = imarissetdataset('uint8', sz(1), sz(2), sz(3), nch, 1);

imarissetvolume(dset, uint8(img(:,:,zsl)),0)

for ch = 2:4

    img2 = imread_bf(filename, struct('series', se, 'time', ti, 'channel', ch, 'x', xr, 'y', yr));
    img2 = imzreverse(squeeze(img2));

    imarissetvolume(dset, uint8(img2(:,:,zsl)),ch-1)
end


%% set colors

vRed = 0.0;
vGreen = 0.0;
vBlue = 1.0;
vAlpha = 0;
vRGBA = [vRed, vGreen, vBlue, vAlpha];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine 
different components (four bytes) into one integer
dset.SetChannelColorRGBA(0, vRGBA);

vRed = 1.0;
vGreen = 1.0;
vBlue = 1.0;
vAlpha = 0;
vRGBA = [vRed, vGreen, vBlue, vAlpha];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine 
different components (four bytes) into one integer
dset.SetChannelColorRGBA(1, vRGBA);

vRed = 1.0;
vGreen = 0.0;
vBlue = 0.0;
vAlpha = 0;
vRGBA = [vRed, vGreen, vBlue, vAlpha];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine 
different components (four bytes) into one integer
dset.SetChannelColorRGBA(2, vRGBA);

vRed = 0.0;
vGreen = 1.0;
vBlue = 0.0;
vAlpha = 0;
vRGBA = [vRed, vGreen, vBlue, vAlpha];
vRGBA = round(vRGBA * 255); % need integer values scaled to range 0-255
vRGBA = uint32(vRGBA * [1; 256; 256*256; 256*256*256]); % combine 
different components (four bytes) into one integer
dset.SetChannelColorRGBA(3, vRGBA);


%% threshold in sox2 / snail intensity

figure(120)
subplot(1,2,1)
hist([stats{3}.MedianIntensity], 255)
xlabel('Snail')

subplot(1,2,2)
hist([stats{4}.MedianIntensity], 255)
xlabel('Sox2')


%%

figure(120)
subplot(1,2,1)
hist([stats{3}.MeanIntensity], 255)
xlabel('Snail')

subplot(1,2,2)
hist([stats{4}.MeanIntensity], 255)
xlabel('Sox2')



%%
%nset = 3;
nset = size(surf);
sfset = surf(1:nset);
fcset = fac(1:nset);
nmset = norm(1:nset);

imarissetsurface('Segmentation', sfset, fcset, nmset);


%%

thsnail = 8.61;
thsox2 = 23.3;


snail = [stats{3}.MeanIntensity];
sox2   = [stats{4}.MeanIntensity];
ncells = length(snail);
snailsox2 = zeros(1, ncells);
snailsox2(sox2 > thsox2) = 2;
snailsox2(snail > thsnail) = 3;
snailsox2(and(snail > thsnail, sox2 > thsox2)) = 1;

imarissetstatistics('snailsox', snailsox2);


%%

% create separate surfaces for different

ids = find(and(sox2 >= thsox2, snail < thsnail));

nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('Sox2', sfset, fcset, nmset);


%%

ids = find(and(snail >= thsnail, sox2 < thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('Snail', sfset, fcset, nmset);


%%
ids = find(and(snail >= thsnail, sox2 >= thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('SnailSox2', sfset, fcset, nmset);

%%
ids = find(and(snail < thsnail, sox2 < thsox2));
nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('None', sfset, fcset, nmset);



%end



%%

save('Z:\aryeh_sox2_snail.mat')