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
fns =  '/var/run/media/ckirst/ChristophsData/Eric/20130911T114815/MIP/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
flds = 1:16;

k = 1;
clear imgs
for f = flds
   img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 1, 'F', f, 'Z', 0)));
   imgs{k} = imfrmtPermute(img, 'XY', 'Yx');
   k = k + 1;
end
imgs = imfrmtReshape(imgs', 'S', 'Uv', 'S', 'UV', [4,4]);
ids = num2cell(1:length(flds), 1)
ids = imfrmtReshape(ids', 'S', 'Uv', 'S', 'UV', [4,4])


%%

figure(1); clf
img = stitchPreview(imgs, 'overlap', 80, 'scale', 0.5);
implot(img)


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

%%

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
fprintf('Alignment: found %g connected components\n', nsubalgn);

if verbose  
   var2char({subalgn.anodes})
end


%%
for s = 1:nsubalgn
   fprintf('\n\nAligning component: %g / %g\n', s, nsubalgn)
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 120, 'overlap.min', 90, 'shift.max', 10);
   if verbose && s < 20 %%&& subalgn(s).nNodes < 75
      subalgn(s).printInfo 
      figure(100+s); clf
      
      subalgn(s).plotPreviewStiched('scale', 0.5)
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



%% Colony Detection 

% detect by resampling and closing
scale = algnAll.source.previewScale;
roi = detectROIsByClosing(algnAll, 'scale', scale, 'threshold', 1.25* th, 'strel', 1, 'radius', 20, 'dilate', 65, 'check', true)


%% Colony 

colonies = Colony(algnAll, roi);
ncolonies = length(colonies);

figure(1); clf
colonies.plotPreview


%% bounding boxes

clc
coloniesC = colonies.copy;
for i = 1:length(colonies)
   roi = colonies(i).roi;

   %[c, r]= fitCircle(roi.p);
   %roic = ROIDisk(c,r);
   roic = roi.boundingBox;

   coloniesC(i).roi = roic;
end

figure(2); clf
coloniesC.plotPreview

%% take full colonies only

figure(8); clf
hist(coloniesC.volume, 20)

coloniesD = coloniesC.copy;
ids = coloniesC.volume > 400000
coloniesD = coloniesD(ids)

% remove lower border
rois = [coloniesD.roi];
pp = [rois.p1];
ids = pp(2,:) > 0;

coloniesD =  coloniesD(ids);


figure(9);
coloniesD.plotPreview


%% Estimate Flatfield / Background

f = 7;
clear imgs
imgf = zeros([512 512]);
imgb = inf * ones([512 512]);

k = 1;

for t = 1:155 %155
   fprintf('t = %d\n', t);
   for f = flds
      img = double(imread(tagExpressionToString(fns, 'C', 2, 'T', t, 'F', f, 'Z', 0)));
      img = imfrmtPermute(img, 'XY', 'Yx');

      imgb = min(imgb, img);
      imgf = ((k-1) * imgf + img)/k;
      
      k = k + 1;
   end
end


%%
figure(13); clf;
implottiling({imgf; imgb});

figure(14); clf;

img = double(imread(tagExpressionToString(fns, 'C', 1, 'T', 16, 'F', 6, 'Z', 0)));
img = imfrmtPermute(img, 'XY', 'Yx');
implot(img)


%% Save

% make sure to clear cache before saving
%coloniesD.clearCache();
%save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/SMAD4_colonies/Analysis/coloniesMIP.mat', 'colonies', 'coloniesD')
%ismip.clearCache();
%save('/home/ckirst/Data/Science/Projects/StemCells/Experiment/SMAD4_colonies/Analysis/sourceMIP.mat', 'ismip')

%% Load
% 
% clear all %#ok<CLSCR>
% close all
% clc
% 
% verbose = true;

%%
%load('/home/ckirst/Data/Science/Projects/StemCells/Experiment/H2BCitrineTrackingAndFate/Analysis/coloniesMIP.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data 


cid = 1;

figure(9);
coloniesD(cid).plotPreview


%%

dd = coloniesD(cid).data;
figure(6); clf;
implot(dd)


%%




%% manual extraction as BF is slow for 1000000 files 

clc

fns =  '/var/run/media/ckirst/ChristophsData/Eric/20130911T114815/Image/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';
texpout = '/data/Science/Projects/StemCells/Experiment/FateDynamics/Sox17/ColonyDataTempAlignNew/ColonyW<Q,2>T<T,4>Z<Z,2>C<C,1>.tif';

cid = 1;

cc = 1;

as = coloniesD(cid).source;
roi = coloniesD(cid).roi;
[ids, shids] = as.nodesFromROI(roi);

shifts = as.imageShifts;
shifts = shifts(shids);

verbose = true;

roiS = roi.copy;
roiS.shift(-shifts{1}+1);
    


for tt = 1:155
   
   fprintf('Time: %d\n', tt)
   
   for cc = 1:2

   for zz = [7, 1:6, 8:15]

      imgs = cell(16,1);
      k = 1;
      for f = flds
         img = double(imread(tagExpressionToString(fns, 'C', cc, 'T', tt, 'F', f, 'Z', zz)));
         imgs{k} = imfrmtPermute(img, 'XY', 'Yx');
         k = k + 1;
      end
      imgs = imfrmtReshape(imgs, 'S', 'Uv', 'S', 'UV', [4,4]);
      %ids = num2cell(1:length(flds), 1)
      %ids = imfrmtReshape(ids', 'S', 'Uv', 'S', 'UV', [4,4])

  
      %% stich and temporally algin

      st = stitchImages(imgs, shifts, 'method', 'Mean'); 
      
      if tt  > 1 && zz == 7 && cc == 1 % should be first zz loop id
         % align with previous image
     
         [sh, q] = align2ImagesByRMS(stpre, st, 'overlap.min', 1200);
      
         fprintf('Time: %d Drift is: %s\n', tt, var2char(sh))
      
         if verbose
            %figure(11); clf
            %implottiling({st; stpre});
            
            figure(22); clf
            plotAlignedImages({stpre; st}, {[0,0]; sh});
            drawnow
         end

         %shifts = cellfunc(@(x) x - sh, shifts);
         %st = stitchImages(imgs, shifts, 'method', 'Mean');
         
         stpre = st;
         
         roiS.shift(-sh);
         
      elseif zz == 7 && cc == 1
         stpre = st;
         sh = [0,0];
      end
        

      
  
      dat = roiS.extractData(st);
      
      if verbose && zz == 7

         figure(20);
         implot(st)
         
         if tt > 1
            figure(13); 
            imsubplot(5,6,tt)
            plotAlignedImages({dat; datpre}, {[0,0]; [0,0]});
            title(num2str(tt));
            drawnow
            datpre = dat;
         else         
            datpre = dat;
         end
            
            
      %size(st)
%          figure(21);
%          implot(dat)
%          drawnow
         
         %max(dat(:))
         %class(dat)
      end
  
      
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
         



