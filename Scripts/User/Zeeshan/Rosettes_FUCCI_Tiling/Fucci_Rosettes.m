%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Citrine Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;

initializeParallelProcessing(12) % number of processors

%%

fns =  '/home/ckirst/Data/Science/Projects/StemCells/Experiment/Fucci/Rosettes/20150323T225416/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
flds = 1:(4*6);
fshp = [6,4];

%flds = 25:44;
%fshp = [5,4];

flds = 45:68;
fshp = [6,4];

flds = 69:88;
fshp = [5,4];

k = 1;
clear imgs
for f = flds
   img = double(imread(tagExpressionToString(fns, 'C', 3, 'T', 1, 'F', f, 'Z', 0)));
   imgs{k} = imfrmtPermute(img, 'XY', 'Yx');
   k = k + 1;
end
imgs = imfrmtReshape(imgs', 'S', 'Uv', 'S', 'UV', fshp);
ids = num2cell(1:length(flds), 1)
ids = imfrmtReshape(ids', 'S', 'Uv', 'S', 'UV', fshp)

%%

figure(1); clf
implottiling(imgs)

%%

figure(1); clf
img = stitchPreview(imgs, 'overlap', 10, 'scale', 0.5);
implot(img)%coloniesD.clearCache();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 

%%
is = ImageSource();
is.fromCell(imgs);

is.printInfo

%%

algn = Alignment(is);
algn.setCaching(true);
algn.printInfo


%% Background Intensity for Overlap Quality of Tiles

d = imgs{2};
figure(2); clf;
imsubplot(2,1,1)
implot(d)


dd = d(:);
size(dd)
dd(dd < 2200) = [];
size(dd)
subplot(1,2,2)
hist(dd(:), 17550)

th = thresholdFirstMin(dd, 'nbins', 250, 'delta', 100)

%%
th = 3500;

figure(30); clf
d =imgs{3};
implot(double(d > th))


%% Quality of Overlap between neighbouring Tiles 

% parameter: see overlapQuality
algn.calculateOverlapQuality('threshold.max', th, 'overlap.max', 100);

hist(algn.overlapQuality, 256)


%% Align components
clc
clear subalgn
subalgn = algn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn);
fprintf('Alignment: foutexpoutnd %g connected components\n', nsubalgn);

if verbose  
   var2char({subalgn.anodes})
end


%%
for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 30, 'overlap.min', 0, 'shift.max', 3);
   if verbose && s < 20 %%&& subalgn(s).nNodes < 75
      subalgn(s).printInfo 
      figure(100+s); clf
      
      subalgn(s).plotPreviewStiched('scale', 1)
   end
end

%% Align Components
clc
subalgn.alignPositions;

% merge to single alignment

algnAll = subalgn.merge;
algnAll.printInfo

%%
if verbose
   figure(5); clf
   algnAll.plotPreviewStiched('scale', 0.5)
end





%% Estimate Flatfield / Background


imgfCh = cell(1,3);
imgbCh = cell(1,3);

parfor ch = 1:3
   imgf = zeros([512 512]);
   imgb = inf * ones([512 512]);
   k = 1;
   for t = 1:180 %155
      fprintf('t = %d\n', t);
      for f = 1:88%flds
         img = double(imread(tagExpressionToString(fns, 'C', ch, 'T', t, 'F', f, 'Z', 0)));
         img = imfrmtPermute(img, 'XY', 'Yx');
         
         imgb = min(imgb, img);
         imgf = ((k-1) * imgf + img)/k;
         
         k = k + 1;
      end
   end
   
   imgfCh{ch} = imgf;
   imgbCh{ch} = imgb;
end


%%

figure(13); clf;
implottiling([cellfunc(@mat2gray, imgfCh); cellfunc(@mat2gray, imgbCh)]);


figure(14); clf;

img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 6, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implot(img)



%%
imgfChs = cellfunc(@(x) filterGaussian(x, 100), imgfCh);
imgbChs = cellfunc(@(x) filterGaussian(x, 100), imgbCh);


figure(13); clf;
implottiling(cellfunc(@mat2gray, imgfChs));


ch = 1;
figure(14); clf;
img = double(imread(tagExpressionToString(fns, 'C', ch, 'T', 16, 'F', 4, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implottiling({mat2gray(img); mat2gray((img - imgbCh{ch})./ (imgfChs{ch} - imgbCh{ch}))})


%%

ch = 3;
imgsc = cell(size(imgs));
parfor j = 1:numel(imgs)
   imgsc{j} = (imgs{j} -imgbCh{ch}) ./ (imgfChs{ch} - imgbCh{ch});
end

is.fromCell(imgsc);

if verbose
   figure(5); clf
   
   algnAll.plotPreviewStiched('scale', 0.4)
end



%% Save

% make sure to clear cache before saving
fnsC =  '/home/ckirst/Data/Science/Projects/StemCells/Experiment/Fucci/Rosettes/20150323T225416_Correction.mat';
save(fnsC, 'imgbCh', 'imgfCh', 'imgfChs')


%% Load
% 

fnsC =  '/home/ckirst/Data/Science/Projects/StemCells/Experiment/Fucci/Rosettes/20150323T225416_Correction.mat';
load(fnsC)

imgfChs = cellfunc(@(x) filterGaussian(x, 100), imgfCh);
imgbChs = cellfunc(@(x) filterGaussian(x, 100), imgbCh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Corrected Movie
mfn = '/home/ckirst/Data/Science/Projects/StemCells/Experiment/Fucci/Rosettes/Area_4.avi';

writerObj = VideoWriter(mfn);
open(writerObj);


ts = [1:183, 185:288]; % 184 is missing ???
%ts = 1:2;

for t = ts
   
   fprintf('movie: %s, time %d / %d\n', mfn, t, 288)
   
   imgC = zeros([algnAll.dataSize, 3]);
   %imgC = zeros([3061, 2049, 3]);
   
   for ch = 1:3
      clear imgs
      k = 1;
      for f = flds
         img = double(imread(tagExpressionToString(fns, 'C', ch, 'T', t, 'F', f, 'Z', 0)));
         imgs{k} = imfrmtPermute(img, 'XY', 'Yx');
         k = k + 1;
      end
      imgs = imfrmtReshape(imgs', 'S', 'Uv', 'S', 'UV', fshp);
      
      imgsc = cell(size(imgs));
      parfor j = 1:numel(imgs)
         imgsc{j} = (imgs{j} - imgbChs{ch}) ./ (imgfChs{ch} -  imgbChs{ch});
      end
      
      is.fromCell(imgsc);
      
      imgC(:,:,ch) = algnAll.stitch('method', 'Max');
   end
   
%    figure(4); clf;
%    implot(imgC)
   
   % rescale
   %
   figure(7); clf
   for c = 1:3
   subplot(1,3,c)
   dd = imgC(:,:,c);
   hist(dd(:), 256)
   end
   
   
   imgsCs = imgC;
   cmin = [0.25, 0.25, 0.3];
   cmax = [9, 14, 12];
   
   for c = 1:3
      imgsCs(:,:,c) = (imclip(imgsCs(:,:,c), cmin(c), cmax(c)) - cmin(c)) / (cmax(c) - cmin(c));
      imgsCs(:,:,c) = imgsCs(:,:,c) - imopen(imgsCs(:,:,c), strel('disk', 20));
   end
   
   %imgsCs = filterAutoContrast(imgsCs);
   
   
   
   figure(8); clf
   implot(imgsCs)
   drawnow
   % write image to movie
   writeVideo(writerObj, imfrmtReformat(imgsCs, 'XYC', 'YXC'));
end

writerObj.close;








%% manual extraction as BF is slow for 1000000 files 

clc

fns =  '/data/Science/Projects/StemCells/Experiment/FateDynamics/Sox17_10min/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '/data/Science/Projects/StemCells/Experiment/FateDynamics/Sox17_10min_Colony/Colony<Q,2>T<T,4>Z<Z,2>C<C,1>.tif';

cid = 2;

cc = 1;

as = coloniesD(cid).sou%%
figure(13); clf;

imgfs = filterGaussian(imgf, 100);

implottiling({mat2gray(imgf);mat2gray(imgb); mat2gray(filterGaussian(imgb, 15)); mat2gray(imgfs)});

figure(14); clf;

img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 10, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implottiling({mat2gray(img); mat2gray((img - imgb)./ (imgfs - imgb))});
roi = coloniesD(cid).roi;
[iids, shids] = as.nodesFromROI(roi)

shifts = as.imageShifts;
shifts = shifts(shids);

verbose = true;
if verbose
   figure(5); clf
   algnAll.plotPreviewStiched('scale', 0.5)
end

roiS = roi.copy;
roiS.shift(-shifts{1}+1);

idsm = cell2mat(ids);

[~, iid] = ismember(iids, idsm(:));
%%
figure(13); clf;

imgfs = filterGaussian(imgf, 100);

implottiling({mat2gray(imgf);mat2gray(imgb); mat2gray(filterGaussian(imgb, 15)); mat2gray(imgfs)});

figure(14); clf;

img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 10, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implottiling({mat2gray(img); mat2gray((img - imgb)./ (imgfs - imgb))})
for tt = 1:72
   
   fprintf('Time: %d\n', tt)
   
   for cc = 1:1

   for zz = 1:15

      imgs = cell(length(iids),1);
      k = 1;
      for f = flds(iid)
         img = double(imread(tagExpressionToString(fns, 'C', cc, 'T', tt, 'F', f, 'Z', zz)));
         imgs{k} = imfrmtPermute(img, 'XY', 'Yx');
         k = k + 1;
      end
      %imgs = imfrmtReshape(imgs, 'S', 'Uv', 'S', 'UV', [4,4]);
      %ids = num2cell(1:length(flds), 1)
      %ids = imfrmtReshape(ids', 'S', 'Uv', 'S', 'UV', [4,4])

  
      %% stich and temporally algin

      st = stitchImages(imgs, shifts, 'method', 'Interpolate'); 
     
      
%       if tt  > 1 && zz == 7 && cc == 1 % should be first zz loop id
%          % align with previous image
%      
%          [sh, q] = align2ImagesByRMS(stpre, st, 'overlap.min', 1200);
%       
%          fprintf('Time: %d Drift is: %s\n', tt, var2char(sh))
%       
%          if verbose
%             %figure(11); clf
%             %implottiling({st; stpre});
%             
%             figure(22); clf
%             plotAlignedImages({stpre; st}, {[0,0]; sh});
%             drawnow
%          end
% 
%          %shifts = cellfunc(@(x) x - sh, shifts);
%          %st = stitchImages(imgs, shifts, 'method', 'Mean');
%          
%          stpre = st;
%          
%          roiS.shift(-sh);
%          
%       elseif zz == 7 && cc == 1
%          stpre = st;
%          sh = [0,0];
%       end
%         

      
  
      
   dat = roiS.extractData(st);
   
       %figure(1); clf
       %implot(dat)
   
   %dat = st;
      
%       if verbose && zz == 7
% 
%          figure(20);
%          implot(st)
%          
%          if tt > 1
%             figure(13); 
%             imsubplot(5,6,tt)
%             plotAlignedImages({dat; datpre}, {[0,0]; [0,0]});
%             title(num2str(tt));
%             drawnow
%             datpre = dat;
%          else         
%             datpre = dat;
%          end
            
            
      %size(st)
%          figure(21);
%          implot(dat)
%          drawnow
         
         %max(dat(:))
         %class(dat)
%      end
  
      
      %% save

      fnout =  tagExpressionToString(texpout, 'Q', cid, 'C', cc, 'T', tt, 'Z', zz);
      imwriteTIFF(int32(dat), fnout)
      
      
   end
      
   end
   
end
      
      %%
%       
%       tdat = imread(fnout);
%       size(tdat)
%       max(tdat(:))
%       min(tdat(:))
% 
%       figure(1); clf
%       implot(tdat)
         



