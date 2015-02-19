%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Citrine Movie - Freds data for Eric / Colony
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

initializeParallelProcessing(12) % number of processors

verbose = true;

%%

texpmip = '/home/siggia/TGMM/Smad4-GFP_trial_quantification/Colony1/W1F<F,3>T<T,4>Z<Z,2>C<C,1>.tif';

ismip = ImageSourceFiles(texpmip);
ismip.printInfo

%%
ismip.printInfo

%%
ismip.resetRange;
ismip.setReshape('F', 'UV', [4,4]);
ismip.setCellFormat('UvTZC');
ismip.printInfo

%%

ismip.setRange('T', 1,  'C', 2)
ismip.printInfo

% %%
% 
% ismip.setRange('T', 1, 'U', 1:4, 'V', 1:3)
% ismip.printInfo
% 
figure(1); clf;
cd = ismip.cell
implottiling(cd)



%%

figure(1); clf;
ismip.plotPreviewStiched('overlap', 80, 'scale', 0.5, 'method', 'Interpolate')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment

% create 
clc
algn = Alignment(ismip, 'UV');
algn.setCaching(false);
algn.printInfo


%% Background Intensity for Overlap Quality of Tiles

d = ismip.data(7);
figure(2); clf;
imsubplot(2,1,1)
implot(d)


dd = d(:);
size(dd)
dd(dd < 2200) = [];
size(dd)
subplot(1,2,2)
hist(dd(:), 150)

th = thresholdFirstMin(dd, 'nbins', 50, 'delta', 100)

%%
th = 2100;

figure(30); clf
d = ismip.data(7);
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
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 95, 'overlap.min', 70, 'shift.max', 10);
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




%% Estimate Flatfield / Background

is = ismip.copy;
is.resetRange;
is.printInfo

%%

c = 1;

fns = texpmip;

imgf = zeros([512 512]);
imgb = inf * ones([512 512]);

k = 1;

times = [12, 18, 20, 25, 30, 35, 40, 60, 90, 120, 150] - 11;
slices = [1];
imgbZ = cell(1,length(slices));
imgfZ = cell(1,length(slices));

for z = slices

   imgf = zeros([512 512]);
   imgb = inf * ones([512 512]);

   k = 1;

   for t = times %119
      fprintf('t = %d\n', t);
      
      imgs = is.cell('T', t, 'Z', z, 'C', c)

      imga = cat(3, imgs{:});
      imgmean = mean(imga, 3);

      imgb = min(imgb, min(imga,[], 3));
      imgf = ((k-1) * imgf + imgmean)/k;
      
      k = k + 1;

   end
   
   imgbZ{z} = imgb;
   imgfZ{z} = imgf;
      
end

%%
clc
for z = slices

   imgfZs{z} = filterGaussian(imgfZ{z}, 20);
   imgbZs{z} = filterGaussian(imgbZ{z}, 20);

   imgf = imgfZ{z};
   imgb = imgbZ{z};
   
   img = is.data('C', 1, 'T', 1, 'Z', z, 'U', 2, 'V', 2);
   
   figure(13 + z); clf;
   implottiling({mat2gray(imgf), mat2gray(imgb); mat2gray(imgfZs{z}), mat2gray(imgbZs{z});
                 mat2gray(img), mat2gray((img - imgb)./ (imgfZs{z} - imgb)); 
                 mat2gray(img -imgbZs{z}) , mat2gray((img - imgbZs{z})./ (imgfZs{z} - imgbZs{z}))
                 }')

end

%%
save('~/Desktop/Smad4Colonies/correction.mat', 'imgfZs', 'imgbZs')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract Data 

texpout = '~/Desktop/Smad4Colonies/Colony1/T<T,4>C<C,1>.tif'

times = [12, 18, 20, 25, 30, 35, 40, 60, 90, 120, 150] - 11;
%times = 1;
slices = [1];
sh = algnAll.imageShifts;

verbose = false;
%verbose= true;

cc = 2;
zz = slices;

parfor t = 1:length(times)
   
   tt = times(t);

   fprintf('Time: %d\n', t)
   
   imgs = is.cell('T', tt, 'C', cc, 'Z', zz);

   %imgsc = imgs;
   %imgsc = cellfunc(@(x) imclip(x, 2010, []), imgs);
   imgsc = cellfunc(@(x) (x - imgbZs{1})./ (imgfZs{1} - imgbZs{1}), imgs);
   
   img = stitchImagesByMean(imgsc, sh, 'method', 'Interpolate');

   if verbose
      figure(1); clf
      hist(img(:), 256)
   
      figure(1); clf
      implot(mat2gray(imclip(img, 0, [])))
   end

   fnout = tagExpressionToString(texpout, 'C', cc, 'T', tt + 11, 'Z', zz-1);
   imwriteTIFF(int32(img / 11 * 2^16), fnout)
end


