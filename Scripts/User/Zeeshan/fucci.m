%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUCCI in Neuronal Rosettes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize()
ilinitialize()


%% Load Data

exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'FUCCI ',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Other/Zeeshan/', ...
                 'ImageDirectoryName', 'FUCCI',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'FUCCI-triple_<channel,s> 25-_s6_t<time>.TIF');

exp.info()

tags = exp.fileTagRange
times = tags.time;
channelnames = tags.channel;

verbose = true;

%% Generate Single Frame


t = 1;


%%

t

imgB = imread(exp.fileName('time', t, 'channel', channelnames{1}))';
imgR = imread(exp.fileName('time', t, 'channel', channelnames{2}))';
imgG = imread(exp.fileName('time', t, 'channel', channelnames{3}))';

imgB = mat2gray(imgB);
imgR = mat2gray(imgR);
imgG = mat2gray(imgG);


img = cat(3, imgR, imgG, imgB);

figure(1); clf;
implot(img)

t = t+1;

%% Load full Movie


timax = 68;
tioff = 6;

imgall = zeros(1344, 1024, timax, 3);

ti = 1;
for t = times((1:3:3*timax)+tioff)
   fprintf('loading image t = %g (%g / %g)\n', t, ti, timax);

   imgB = imread(exp.fileName('time', t, 'channel', channelnames{1}))';
   imgR = imread(exp.fileName('time', t, 'channel', channelnames{2}))';
   imgG = imread(exp.fileName('time', t, 'channel', channelnames{3}))';

   %imgB = mat2gray(imgB);
   %imgR = mat2gray(imgR);
   %imgG = mat2gray(imgG);

   imgall(:,:,ti,1) = imgR; imgall(:,:,ti,2) = imgG; imgall(:,:,ti,3) = imgB;

   ti = ti + 1;
end
   
% scale

for c = 1:3
   imgall(:,:,:,c) = mat2gray(imgall(:,:,:,c));
end


%% scale red 

for c = 1
   imgall(:,:,:,c) = 2 * imgall(:,:,:,c);
end

% %%
% for c = 3
%    imgall(:,:,:,c) = 0 * imgall(:,:,:,c);
% end
% 
% %%
% for c = 1:2
%    imgall(:,:,:,c) = 0 * imgall(:,:,:,c);
% end


%%

figure(1); clf;

for ti = 1:timax
   implot(squeeze(imgall(:,:,ti,:)))
   drawnow
   pause(0.4)
end

%%

% try to detect cells

% use nuclear + green + red
ti = 1;
img = 1 * imgall(:,:,ti,1) + imgall(:,:,ti,2) + 1 * imgall(:,:,ti,3);
img = mat2gray(img);

figure(1); clf; imcolormap('b')
implot(img)

%% mask

imgmask = img > 0.15;

figure(1); clf;
implot(imoverlaylabel(img, double(imgmask)))

%% Seeds

imgf =img;
imgf = logFilter(max(imgf(:)) - imgf, [25,25], []);
imgf = mat2gray(imgf);


% h-max detection (only local maxima with height > hmax are considered as maxima
imgmax = imextendedmax(mat2gray(imgf),  0.005);

%local max
%imgmax = imregionalmax(imgf);

% constrain to maxima within mask
imgmax = immask(imgmax, imgmask);

% Combine nearby points
imgmax = imdilate(imgmax, strel('disk', 2));
% fill holes - combination of nearby points can lead to holes
imgmax = imfill(imgmax,'holes');
% shrink to single points - extended maxima usually give better segmentation results
% imgmax = bwmorph(imgmax,'shrink',inf);           

% plot the results.
if verbose  
   figure(20)
   set(gcf, 'Name', ['Seeding'])
   implottiling({imoverlay(img, imgmax), imoverlay(imgf, imgmax)});
end


%% Label

imglab = bwlabeln(imgmax);

figure(21)
implot(imoverlaylabel(img, imglab))


%% Classify Labels

stats = imstatistics(imglab, 'PixelIdxList');


statsR = imstatistics(imglab, stats, 'MeanIntensity', imgall(:,:,ti, 1));
statsG = imstatistics(imglab, stats, 'MeanIntensity', imgall(:,:,ti, 2));

[~, class] = max(cat(2, 1.5 * [statsR.MeanIntensity]', [statsG.MeanIntensity]'), [], 2);

if verbose 
   figure(99); clf
   subplot(1,3,1)
   hist([statsR.MedianIntensity])
   subplot(1,3,2)
   hist([statsG.MedianIntensity])
   subplot(1,3,3)
   hist(class)

   % plot classification 

   pixl = {stats.PixelIdxList};

   imgcls = imglab;
   for i = 1:length(pixl)
      imgcls(pixl{i}) = class(i);
   end

   figure(100); clf;
   imgp = imoverlay(squeeze(imgall(:,:,ti,:)), imgcls == 1, 'r');
   imgp = imoverlay(imgp, imgcls == 2, 'g');
   
   implot(imgp);
end








