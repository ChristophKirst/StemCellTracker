%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analyze Ca Dynamics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Correlation in Fast dynamics (w.r.t some referrnce region)

rpath =  '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/10_54_fast_dynamics/';

%fh = FileHandler('BaseDirectory', rpath, ...
%                  'ImageDirectoryName', '.',...
%                  'ReadImageCommandFormat', 'imread(''<file>'')',...
%                  'ReadImageFileFormat', 'laser40_EM1000_exp100_30frames_001_t<time, d>.TIF');
               
               
rpath =  '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/11_05_fast_dynamics/';

fh = FileHandler('BaseDirectory', rpath, ...
                  'ImageDirectoryName', '.',...
                  'ReadImageCommandFormat', 'imread(''<file>'')',...
                  'ReadImageFileFormat', 'laser50_EM1000_exp100_60frames_001_t<time, d>.TIF');

               

fh.info()
iif = imfinfo(fh.fileName('time', 1));
isize = [iif.Width, iif.Height];
trange = fh.fileTagRange;
trange = trange.time;
datsize = [isize, max(trange)]

%%

mmaquisitiontime(fh.fileName('time', 1))
mmaquisitiontime(fh.fileName('time', 2))


%%

dat = zeros(datsize);
for t = trange
   dat(:,:,t) = fh.readImage('time', t);
end
dat = double(dat);
size(dat)

%figure(1); clf
%implottiling(dat)

%%
nbin = 4; %% pixel bin
img = imbin(dat, [nbin, nbin, 1]);

size(img)

figure(2); clf;
implottiling(img)

%%

figure(3); clf; hold on
cmap = jet(numel(img(:,:,1)));
sp = 2;
for i = 1:sp:size(img, 1)
   for j = 1:sp:size(img, 2)
      plot(diff(squeeze(img(i,j,:))), 'Color', cmap(i+(j-1)*size(img,1),:))
   end
end

%%

figure(3); clf; hold on
cmap = jet(numel(img(:,:,1)));
for i = 1:1
   for j = 1:1
      plot(squeeze(img(i,j,:)), 'Color', cmap(i+(j-1)*size(img,1),:))
   end
end




%% time scales


rpath = '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_17/GCamp6s_Seeding_Rock/';

fn = dir(fullfile(rpath, '*.tif'))


imfinfo(fullfile(rpath, fn(1).name))


%% time scales

rpath = '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_03_27/01_Ca/p0/';

fn = dir(fullfile(rpath, '*.tif'))

iff = imfinfo(fullfile(rpath, fn(1).name))
