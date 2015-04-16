%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brain Slices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

initialize()
%%

bfinitialize
initializeParallelProcessing(12)

addpath('./Scripts/User/Zeeshan/BrainSlices/');

datadir = '/data/Science/Projects/StemCells/Experiment/BrainSections/';
dataname = 'PCW10 Sample A Slide 53 SATB2 TBR1 BRN2 CTIP2';

verbose = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Source %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%is = ImageSourceBF('/data/Server/smb/upload/Brain sections/Sample B Slide 33 SATB2 CUX2 NURR1 CTIP2.lsm');
is = ImageSourceBF([datadir  dataname '.czi']);
clc
is.printInfo

%%

is.setReshape('S', 'UV', [2, 2]);
is.setCellFormat('Uv')
is.setRange('C', 1);
clc; is.printInfo

%%

imgs =  is.cell('Z', 5, 'C', 1);
figure(8); 
implottiling(imgs)

%%

figure(1); clf
is.setRange('C', 1, 'Z', 5);
is.plotPreviewStiched('overlap', 50, 'scale', 0.1);


%% DAPI ZO1

zs = [4,6,8,10,11];
xy = {225:425, 250:450};
si = [cellfun(@length, xy), 3];

ch = [1,2];
rgb = [3,2];
cl = {[2000, 4095], [478, 4095]};

clear imgs
for z = 1:length(zs)
   img = zeros(si);
   for c = [1,2]
      imgC = is.data('Z', zs(z), 'C', ch(c), 'U', 1, 'V', 1);
      imgC = imclip(imgC, cl{c}(1), cl{c}(2));
      imgC = (imgC - cl{c}(1)) /(cl{c}(2) - cl{c}(1));
      img(:,:,rgb(c)) = imgC(xy{:});
   end
   imgs{z} = img;
end
imgs{end+1} = max(cat(4, imgs{:}), [], 4);

h = figure(2); clf;
for i = 1:6
   imsubplot(3,2,i)
   implot(imgs{i})
   axis off
   xlabel(''); ylabel('')
   if i < 6
      imcolorlabel({sprintf('z=%0.1f um', (zs(i)*0.87))}, {[1,1,1]}, 'location', 'sw')
   else
      imscalebar('color', [1,1,1], 'unit', 'um', 'location', 'se', 'length', 10, 'scale', 0.83, 'textoffset', 0.03)
      imcolorlabel({'DAPI', 'ZO1'},{[0,0,1], [0,1,0]}, 'location', 'sw')
   end
end


%%
saveas(h, [datadir dataname '_ZO1' '.pdf']);




%% DAPI PKC

zs = [4,6,8,10,11];
xy = {225:425, 250:450};
si = [cellfun(@length, xy), 3];

ch = [1,3];
rgb = [3,1];
cl = {[2000, 4095], [478, 4095]};

clear imgs
for z = 1:length(zs)
   img = zeros(si);
   for c = [1,2]
      imgC = is.data('Z', zs(z), 'C', ch(c), 'U', 1, 'V', 1);
      imgC = imclip(imgC, cl{c}(1), cl{c}(2));
      imgC = (imgC - cl{c}(1)) /(cl{c}(2) - cl{c}(1));
      img(:,:,rgb(c)) = imgC(xy{:});
   end
   imgs{z} = img;
end
imgs{end+1} = max(cat(4, imgs{:}), [], 4);

h = figure(2); clf;
for i = 1:6
   imsubplot(3,2,i)
   implot(imgs{i})
   axis off
   xlabel(''); ylabel('')
   if i < 6
      imcolorlabel({sprintf('z=%0.1f um', (zs(i)*0.87))}, {[1,1,1]}, 'location', 'sw')
   else
      imscalebar('color', [1,1,1], 'unit', 'um', 'location', 'se', 'length', 10, 'scale', 0.83, 'textoffset', 0.03)
      imcolorlabel({'DAPI', 'PKC'},{[0,0,1], [1,0,0]}, 'location', 'sw')
   end
end


%%
saveas(h, [datadir dataname '_PKC' '.pdf']);


