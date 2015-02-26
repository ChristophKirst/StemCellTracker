%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nanong on TGMM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test image for mat needed -> make uint16 tif stack using imwrite(img(:,:,z), 'WriteMode', 'append', 'Compression', 'none');
% 
% img = imread('/home/siggia/TGMM/Sox17/data/T0001/ColonyW01T_T0001_C1.tif');
% size(img)
% 
% class(img)

addpath('/home/siggia/matlab/')

%%
clear all
clear classes
close all
clc
initialize
bfinitialize
ijinitialize

verbose = true;

addpath('./Scripts/User/Christoph/mESCNanog');
addpath('./Interface/User/EmbryoData');


initializeParallelProcessing(12) % number of processors

%plotslices = 1:2:40;
tiling = [8,5];



%% Data

datadir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/12Aug14FGFonNanongH2BGFP-700/';
%dexp = 'T<F,3>/T<F,3>.tif';
%fulldexp = fullfile(datadir, dexp);
fulldexp = fullfile(datadir, '12Aug14FGFonNanogH2BGFP-700_movie.lsm');
dns = tagExpressionToFiles(fulldexp);

isd = ImageSourceBF(fulldexp);
isd.printInfo


%% Plot Data
f = 1;

img = isd.data('T', f, 'C', 1);

figure(1); clf
implottiling(img, 'tiling', tiling)

%% Create Filtered Data in TGMM data structure

basedir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/TGMM/data';
mkdir(basedir)

timedir = fullfile(basedir, 'T<T,4>');
outfile = fullfile(timedir, 'T<T,4>.tif');

tmax = 20;

%%
parfor t = 1:tmax
   
   % load data
   img = isd.data('T', t, 'C', 1); 
   
   %imax = 65;
   imax = 255;
   img = imclip(img, 0, imax);
   img = double(img) / imax;
   %img = imclose(img, strel('disk', 1));
   %img = filterMedian(img, [4,4,2]);
   
%    figure(1); clf; 
%    hist(img(:), 512);
%    
%    figure(2); clf;
%    implottiling(img, 'tiling', tiling)
   
   
   %%
   %img = imclip(img, 0, 1);
   %imgf = filterBM(img, 'sigma', 5, 'depth', 3);  % sigma is noise level out of 0..255, %img should be 0..1 however
   %imgf = double(imgf);
   imgf = img;
   
   %figure(2); clf;
   %implottiling(imgf, 'tiling', tiling)
   
   
%    
% 
%    figure(3); clf;
%    imgp = zeros([size(img,1), size(img,2), 2* size(img,3)]);
%    k = 1;
%    for i = 1:size(img,3)
%       imgp(:,:,k) = img(:,:,i);
%       imgp(:,:,k+1) = imgf(:,:,i);
%       k = k + 2;
%    end
%    implottiling(imgp, 'tiling', [16,5]);
    
   %% write to stack
   mkdir(tagExpressionToString(timedir, 'T', t-1));
   fn = tagExpressionToString(outfile, 'T', t-1);
   imwriteStack(uint16(imgf * 2^14), fn, 'Compression', 'none');
end


%% Check Data
% 
fn = tagExpressionToString(outfile, 'T', t);
imgr = imreadBF(fn);
max(imgr(:))
min(imgr(:))

figure(5); clf
implottiling(imgr, 'tiling', tiling)


% 
% 
% figure(6); clf
% hist(imgr(:), 256)
% 
% 
% %%
% imgrr = imgr;
% imgrr(imgrr < 200) = 0;
% 
% figure(system(['rm -r /data/Science/Projects/StemCells/Experiment/Mouse/Nanog/TGMM/output/*']);7); clf
% implottiling(imgrr, 'tiling', tiling)



%% Run TGMM 

% in console 
% parallel -j12 /home/ckirst/Programs/TGMM/build/nucleiChSvWshedPBC/ProcessStack TGMM_configFile.txt {1} ::: {0..9} 
% /home/ckirst/Programs/TGMM/build/TGMM  configFile.txt 0 30 > log.txt


%% Clean TGMM

tmax = 20
for t = 1:tmax
   fn = tagExpressionToString(fullfile(timedir, '*.bin'), 'T', t-1)
   delete(fn)
   
   fn = tagExpressionToString(fullfile(timedir, '*.txt'), 'T', t-1)
   delete(fn)
end
   
system(['rm -r /data/Science/Projects/StemCells/Experiment/Mouse/Nanog/TGMM/output/*']);


%% Parsing the Supervoxel Result



%%






