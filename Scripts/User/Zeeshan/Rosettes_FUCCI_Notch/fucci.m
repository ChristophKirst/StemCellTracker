%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUCCI in Neuronal Rosettes + Notch Inhibition %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make tif files consistent
%rename '4CFP 25-' 4CFP *                                                                                                                       
%rename '3RFP 25-' 3RFP *
%rename '2YFP 25-' 2YFP *

clear all
close all

initialize()
bfinitialize()

addpath('./Scripts/User/Zeeshan/Rosettes_FUCCI_Notch');

initializeParallelProcessing(12) % number of processors

conditions = {'Control', 'DAPT'};

expids.(conditions{1}) = [6,8,10];
expids.(conditions{2}) = [3,7,8];

sampleid.(conditions{1}) = '3';
sampleid.(conditions{2}) = '4';


%%
for cond = conditions

fprintf('Running Condition: %s\n', cond{1})
%%

%imgdir = '/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/Control';
%resdir = '/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/Results/Control';

%fnsexp =[imgdir, '/FUCCI_Sample3_w<C,s,4>_s<S>_t<T>.TIF'];

%%
%imgdir = '/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/DAPT';
%resdir = '/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/Results/DAPT';

%fnsexp =[imgdir, '/FUCCI_Sample4_w<C,s,4>_s<S>_t<T>.TIF'];
%tagExpressionToTags(fnsexp)

%%

imgdir = ['/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/', cond{1}];
resdir = ['/data/Science/Projects/StemCells/Experiment/Fucci/FUCCI_DAPT_2014-12-29/Results/'];
resimgdir = [resdir, cond{1}];

fnsexp =[imgdir, '/FUCCI_Sample', sampleid.(cond{1}), '_w<C,s,4>_s<S>_t<T>.TIF'];

%%

is = ImageSourceFiles(fnsexp);
is.printInfo


%%
if false 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Single Frame
t = 1;
pos = 6; % 6, 8, 10 Control
pos = 3; % 3, 7, 8% DAPT

%%
img = fucci_load(is, t, pos);

figure(1); clf;
implot(img)

t = t+1
   img2 = img(:,:,2) - imopen(img(:,:,2), strel('disk', 50));
%%

t = 1;

img = fucci_load(is, t, pos);

imgR = imgray2color(img(:,:,1), 'r');
imgG = imgray2color(img(:,:,2), 'g');
imgB = imgray2color(img(:,:,3), 'b');

figure(3); clf
implottiling({imgB, imgG, imgR;
              imoverlayalpha(imgdic, imgB, 1.2), imoverlayalpha(imgdic, imgG), imoverlayalpha(imgdic, imgR, 1.5)});

           
figure(4); clf
implottiling(imoverlayalpha(imgdic, imgG + 1.2 * imgR, 1.2))


%%
imgovl = imoverlayalpha(imgdic * 1.2 , 0.8 * imgG + 1.25 * imgR, 1.75);
figure(4); clf
implot(imgovl)
imwrite(impqlpermute(imgovl, 'pqc', 'qpc'), fullfile(resdir, ['rosette_overlay_T', num2str0(t, 3), '.tif']))


%%
           
t = t+1;
%% Pseudo Movie

figure(1); clf;

tt = 1:20:217;
timax = length(tt);

imgall = zeros(1343, 1024, timax, 3);

ti = 1;
for t = tt
   imgall(:,:,ti,:) = permute(fucci_load(is,t, pos), [1,2,4,3]);
   implot(squeeze(imgall(:,:,ti,:)))
   drawnow
   ti = ti+1;
end

%%
for ti = 1:timax
   implot(squeeze(imgall(:,:,ti,:)))
   drawnow
   pause(0.4)
end

end


%%
for pos = expids.(cond{1})


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segment Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

tstart = 1;
tend = 217;
tdelta = 1;

texclude = [];

tt = tstart:tdelta:tend;
%tt = [1];

for t = texclude
 k = find(tt == t);
 tt(k) = [];
end

%tt = [tt(1), tt(end)];
%tt = [1];

timax = length(tt)

verbose = false;

%% Loop over times

statssave = repmat({},timax,1);

parfor ti = 1:timax
   t = tt(ti);

   fprintf(['\ncondition: ', cond{1}, ' time: %g  index: %d / %d\n'], t, ti, timax); 
   
   img = fucci_load(is, t, pos);

   % Segmentation and Classification

   % use nuclear + green + red

   
   %imgseg = 1 * img(:,:,1) + 1 * img(:,:,2) + 1 * img(:,:,3);
   img1i = img(:,:,1) - imopen(img(:,:,1), strel('disk', 50));
   img2i = img(:,:,2) - imopen(img(:,:,2), strel('disk', 50));
   img3i = img(:,:,3) - imopen(img(:,:,3), strel('disk', 50));
   imgseg = 1 * img1i + 1 * img2i + 0 * img3i;

   imglab = fucci_label(imgseg, verbose);
   
   [stats, imgl, img3] = fucci_classify(img(:,:,:), imglab, imgseg, verbose);

   %%
   % Movie
   pixl = {stats.PixelIdxList};
   class = [stats.class];

   imgcls = imglab;
   for i = 1:length(pixl)
      imgcls(pixl{i}) = class(i);
   end

   imgp = imoverlay(img3(:,:,:), imgcls == 1, 'r');
   imgp = imoverlay(imgp, imgcls == 2, 'g');
   imgp = imoverlay(imgp, imgcls == 3, 'b');
   
   imgdic = is.data('T', t, 'S', pos, 'C', 1);
   imgs2 = imgray2color(mat2gray(imgdic), 'w');
   imgs2 = imoverlayalpha(imgs2, imgray2color(img3(:,:,2), 'g') + 1.5 * imgray2color(img3(:,:,1), 'r'), 3);

   imgs2 = imoverlay(imgs2(:,:,:), imgcls == 1, 'r');
   imgs2 = imoverlay(imgs2, imgcls == 2, 'g');
   %imgs2 = imoverlay(imgs2, imgcls == 3, 'b');

   %save movie
   imwriteBF(imgs2,  fullfile(resimgdir, [cond{1}, '_classification_S', num2str0(pos,2), '_T', num2str0(ti, 3), '.tif']))

   if verbose
      figure(100); clf; 
      implottiling({imgp; imgs2})

      figure(101); clf;
      hist(class)
   end
   
   %drawnow;
  
   statssave{ti} = stats;
   

   %ti = ti +1;
 
end


%% Save the Results to Disk

save(fullfile(resdir, ['stats_', cond{1}, '_', num2str(pos), '.mat']), 'statssave');


%% Load Data

% load(fullfile(resdir, ['stats_', cond{1}, '_', num2str(pos), '.dat']))


%% Plot Time Evolution

nt =  length(statssave);

cs = zeros(nt, 3);

for t = 1:nt
   cs(t, 1) = total([statssave{t}.class] == 1);
   cs(t, 2) = total([statssave{t}.class] == 2);
   cs(t, 3) = cs(t,1) + cs(t,2);
   
end

darkgreen = [0, 0.5, 0];

h = figure(105); clf; 

subplot(1,2,1)

hold on
plot(cs(:,2), 'color', darkgreen)
plot(cs(:,1), 'r')

subplot(1,2,2); hold on
plot(cs(:,2)./cs(:,3), 'color', darkgreen)
plot(cs(:,1)./cs(:,3), 'color', 'r')

%%
print(h,'-dpdf', fullfile(resdir,['result_', cond{1}, '_', num2str(pos), '.pdf']))

%% save and plot with mathematica

save(fullfile(resdir, ['result_', cond{1}, '_' , num2str(pos), '_mathematica.mat']), 'cs', '-v4')


end % positions

end % conditions 




%% Plot

pn = 1;
for cond = conditions

   for pos =  expids.(cond{1})


%% Load Data

load(fullfile(resdir, ['stats_', cond{1}, '_', num2str(pos), '.mat']))


%% Plot Time Evolution

nt =  length(statssave);

cs = zeros(nt, 3);

for t = 1:nt
   cs(t, 1) = total([statssave{t}.class] == 1);
   cs(t, 2) = total([statssave{t}.class] == 2);
   cs(t, 3) = cs(t,1) + cs(t,2);
   
end

darkgreen = [0, 0.5, 0];

h = figure(105); 
subplot(2,3,pn)
hold on
plot(cs(:,2), 'color', darkgreen)
plot(cs(:,1), 'r')
ylim([0,900]);

h2 = figure(106); 
subplot(2,3,pn); hold on
plot(cs(:,2)./cs(:,3), 'color', darkgreen)
plot(cs(:,1)./cs(:,3), 'color', 'r')
ylim([0,1]);


pn = pn + 1;

   end
end



%%
print(h,'-dpdf', fullfile(resdir,['result_count_all.pdf']))
print(h2,'-dpdf', fullfile(resdir,['result_fraction_all.pdf']))



%% save and plot with mathematica

save(fullfile(resdir, ['result_', cond{1}, '_' , num2str(pos), '_mathematica.mat']), 'cs', '-v4')




