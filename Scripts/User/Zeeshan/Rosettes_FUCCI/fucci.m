%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUCCI in Neuronal Rosettes %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

initialize()

addpath('./Scripts/User/Zeeshan');



%% Load Data

exp = Experiment('name', 'Example', 'description', 'FUCCI ',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Other/Zeeshan/', ...
                 'ImageDirectoryName', 'FUCCI',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'FUCCI-triple_<channel,s>_s6_t<time>.TIF');

exp.info()

tags = exp.findImageTagRange
times = tags.time;
channelnames = tags.channel;

verbose = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Single Frame
t = 1;

%%
img = fucci_load(exp, t);

figure(1); clf;
implot(img)

t = t+1

%%

t = 217

img = fucci_load(exp, t);

imgdic = exp.readImage('time', t, 'channel', 'w1DIC')';
imgdic = mat2gray(imgdic);
imgdic = gray2rgb(imgdic);

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
imwrite(impqlpermute(imgovl, 'pqc', 'qpc'), fullfile(exp.ResultDirectory, ['rosette_overlay_T', num2str0(t, 3), '.tif']))




%%
           
t = t+1;


%% Pseudo Movie

figure(1); clf;

tt = 1:100:217;
timax = length(tt);

imgall = zeros(1344, 1024, timax, 3);

ti = 1;
for t = tt
   imgall(:,:,ti,:) = permute(fucci_load(exp,t), [1,2,4,3]);
   implot(squeeze(imgall(:,:,ti,:)))
   ti = ti+1;
end

%%
for ti = 1:timax
   implot(squeeze(imgall(:,:,ti,:)))
   drawnow
   pause(0.4)
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segment Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

rc = exp.loadData('rosette_centers.mat');

tstart = 1;
tend = 217;
tdelta = 1;

texclude = [];

tt = tstart:tdelta:tend;

for t = texclude
 k = find(tt == t);
 tt(k) = [];
end

%tt = [tt(1), tt(end)];
%tt = [1];

timax = length(tt);
tt

clear fstats
clear result

%% Loop over times
ti = 1;
for t = tt

   fprintf('time: %g  index: %d / %d\n', t, ti, timax); 
   
   img = fucci_load(exp, t);

   % Segmentation and Classification

   % use nuclear + green + red

   
   %imgseg = 1 * img(:,:,1) + 1 * img(:,:,2) + 1 * img(:,:,3);
   imgseg = 0 * img(:,:,1) + 1 * img(:,:,2) + 1 * img(:,:,3);

   imglab = fucci_label(imgseg, verbose);


   [stats, imgp] = fucci_classify(img(:,:,:), imglab, imgseg, verbose);

   % Rosette Center
   rcenter = fix(rc(t, (1:2)));

   % Radial Statistics
   stats = fucci_stats(imglab, rcenter, stats, verbose);

   
   fstats{ti} = stats;

   % Measure
   
   isize = size(imglab);
   rmax = 700; rmin = 100;
   dr = 50;
   rr = rmin:dr:rmax;
   rmax = rmax + dr;
   
   th1 = 3/4*pi - 0.1;
   th2 = 5/4 * pi;

   ri = 1;
   clear cR cG cRr cGr ar
   cR= zeros(1, length(rr));
   cG = cR; cRr = cR; cGr = cR; ar = cR;
   
   allid = [];
   for r = rr
      %count green / red
      
      id = and([stats.dist] >= r, [stats.dist] < r+dr);
      id = and(id, and([stats.theta] >= th1 , [stats.theta] <= th2));
      
      allid = [find(id), allid];
      
      cR(ri) = sum([stats(id).class] == 1);
      cG(ri) = sum([stats(id).class] == 2);
      
      %{sum(id), cR(ri), cG(ri)}
      
      ar(ri) = overlapDiskSegmentRectangle(rcenter, r + dr, th1, th2, [0,0], isize) - overlapDiskSegmentRectangle(rcenter, r, th1, th2, [0,0], isize);
      
      cRr(ri) = cR(ri) / ar(ri);
      cGr(ri) = cG(ri) / ar(ri);
      
      ri = ri +1;
   end


   if verbose
      figure(54); clf;
      
      subplot(2,1,1)
      hold on
      plot(rr, cRr, 'r*-')
      plot(rr, cGr, 'g*-')
      
      subplot(2,1,2)
      hold on
      plot(rr, cR, 'r*-')
      plot(rr, cG, 'g*-')
   end
   
   result(ti).rr = rr;
   result(ti).cR = cR;
   result(ti).cRr = cRr;
   result(ti).cG = cG;
   result(ti).cGr = cGr;  
   result(ti).ar = ar;
   result(ti).id = allid;  
   
 
   % Movie

   pixl = {stats(allid).PixelIdxList};
   class = [stats(allid).class];

   imgcls = imglab;
   for i = 1:length(pixl)
      imgcls(pixl{i}) = class(i);
   end

   imgp = imoverlay(img(:,:,:), imgcls == 1, 'r');
   imgp = imoverlay(imgp, imgcls == 2, 'g');
   imgp = imoverlay(imgp, imgcls == 3, 'b');
   
   imgs2 = imgray2color(mat2gray(imgseg), 'w');
   imgs2 = imoverlay(imgs2(:,:,:), imgcls == 1, 'r');
   imgs2 = imoverlay(imgs2, imgcls == 2, 'g');
   imgs2 = imoverlay(imgs2, imgcls == 3, 'b');

   %add center as white dot   
   dotr = 10;
   
   irc = fix(rcenter);
   [ii, jj] = meshgrid(-dotr:dotr,-dotr:dotr);
   ii = ii(:); jj = jj(:);
   pp = [ii, jj]; %length(pp)
   pp = pp(ii.^2 + jj.^2 < dotr^2, :); %length(pp)
   pp = pp + repmat(irc, length(pp), 1);
   

   imgp2 = imgp;

   
   
   for k = 1:length(pp)
      imgp2(pp(k,1), pp(k,2),:) = 255;
      imgs2(pp(k,1), pp(k,2),:) = 255;
   end
   
   for c = 1:3
      p2 = fix(- rmax * [1, tan(th1)] + rcenter);
      imgp2(:,:,c) = impixelline(imgp2(:,:,c), rcenter, p2, 255);
      imgs2(:,:,c) = impixelline(imgs2(:,:,c), rcenter, p2, 255);
      
      
      p2 = fix(- rmax * [1, tan(th2)] + rcenter);
      imgp2(:,:,c) = impixelline(imgp2(:,:,c), rcenter, p2, 255);
      imgs2(:,:,c) = impixelline(imgs2(:,:,c), rcenter, p2, 255);
      
      imgp2(:,:,c) = impixelcircle(imgp2(:,:,c), rcenter, rmin, th1, th2, 255);
      imgp2(:,:,c) = impixelcircle(imgp2(:,:,c), rcenter, rmax, th1, th2, 255);
      
      imgs2(:,:,c) = impixelcircle(imgs2(:,:,c), rcenter, rmin, th1, th2, 255);
      imgs2(:,:,c) = impixelcircle(imgs2(:,:,c), rcenter, rmax, th1, th2, 255);
   end

   figure(101); clf;
   implottiling({imgp2, imgs2; img, imgray2color(mat2gray(imgseg), 'w')})
   
   imwrite(impqlpermute(imgp2, 'pqc', 'qpc'), fullfile(exp.ResultDirectory, ['rosette_T', num2str0(ti, 3), '.tif']))

   ti = ti +1;
 
end


%% Save the Results to Disk

exp.saveData('result.mat', result);
exp.saveData('stats.mat', fstats);



%% Load Data

result = exp.loadData('result.mat');
fstats = exp.loadData('stats.mat');



%% Dist Cell Center vs. Time and Class

% for each cell plot radius on y and in color the class

ti = 1;
timax  = length(fstats);

rmin = 100;
rmax = 850;

th1 = 3/4*pi - 0.1;
th2 = 5/4 * pi;


tdat = [];
ddat = [];
cdat = [];
rdat = [];
gdat = [];

for ti = 1:1:timax
   fprintf('time %d / %d\n', ti, timax);
     
   stats = fstats{ti};
   
   id = and([stats.dist] >= rmin, [stats.dist] <= rmax);
   id = and(id, and([stats.theta] >= th1 , [stats.theta] <= th2));

   tdat = [tdat, ti* ones(1,sum(id))];
   ddat = [ddat, [stats(id).dist]];
   rdat = [rdat, [stats(id).RIntensity]];
   gdat = [gdat, [stats(id).GIntensity]];
   cdat = [cdat, [stats(id).class]];
   
end

%%
h = figure(17); clf;
hold on

ms = 5;

id = (cdat == 2);
id = and(id, ddat < 850);
plot(tdat(id), ddat(id), 'g.', 'MarkerSize', ms);

id = (cdat == 1);
id = and(id, ddat < 850);
plot(tdat(id), ddat(id), 'r.', 'MarkerSize', ms);

id = (cdat == 3);
id = and(id, ddat < 850);
%plot(tdat(id), rdat(id), 'b.');

xlim([0,215])
ylim([0, 900])



%%
print(h,'-dpdf', exp.ResultFile(['resotte_radius_cells.pdf']))


%%

exp.saveData('mathematica.mat', [tdat; ddat; cdat; gdat; rdat], '-v4')

%% Intensity vs. Time

% for each cell plot radius on y and in color the class

ti = 1;
rmax = 750;


th1 = 3/4*pi - 0.1;
th2 = 5/4 * pi;


rdat = [];
gdat = [];
tdat = [];
ddat = [];
cdat = [];

for ti = 1:timax
   ti

   stats = fstats{ti};

   tdat = [tdat, ti* ones(1,length([stats.dist]))];
   rdat = [rdat, [stats.RIntensity]];
   gdat = [gdat, [stats.GIntensity]];
   ddat = [ddat, [stats.dist]];
   cdat = [cdat, [stats.class]];
end



%%
h = figure(17); clf;
hold on

ids = 1:length(tdat);
ids = 1:100;

scatter(tdat(ids), ddat(ids), 20, squeeze(imgray2color(rdat(ids), 'r')), 'fill');




%%
plot(tdat, rdat, 'r.');
xlim([0,215])

%%
print(h,'-dpdf', exp.ResultFile(['resotte_intensity_cells.pdf']))

%% save and plot with mathematica

exp.saveData('mathematica.mat', [tdat; ddat; dat; gdat], '-v4')




%% Area Normalization and Binning

rc = exp.loadData('rosette_centers.mat');
tmax = size(rc,1)
isize = [1344, 1024];


%% r in um !!!
rmin = 30;
rmax = 240; 
dr = 30;
rr = rmin:dr:rmax;
length(rr)
rmax = rmax + dr;

%%
clear res
for t = 1:tmax
   ti = t;

   % Statistics
   stats = fstats{ti};
   
   % Rosette Center
   rcenter = fix(rc(t, (1:2)));

   th1 = 3/4*pi - 0.1;
   th2 = 5/4 * pi;

   clear cR cG cRr cGr ar
   cR= zeros(1, length(rr));
   cG = cR; cRr = cR; cGr = cR; ar = cR;

   ri = 1;
   for r = rr  
      id = and([stats.dist] >= r / 0.32, [stats.dist] < (r+dr)/0.32);
      id = and(id, and([stats.theta] >= th1 , [stats.theta] <= th2));
      
      cR(ri) = sum([stats(id).class] == 1);
      cG(ri) = sum([stats(id).class] == 2);
   
      %{sum(id), cR(ri), cG(ri)}
      
      ar(ri) = overlapDiskSegmentRectangle(rcenter, (r + dr)/0.32, th1, th2, [0,0], isize) - overlapDiskSegmentRectangle(rcenter, r/0.32, th1, th2, [0,0], isize);
      
      cRr(ri) = cR(ri) / ar(ri);
      cGr(ri) = cG(ri) / ar(ri);
      
      ri = ri +1; 

   end
   
   %nr(t, :) = ar'; 
   res(t, :) = [ar, cR, cG, cRr, cGr];
   
end


%%

exp.saveData('rosette_statistics.mat', res, '-v4')











%% Class as Density Plots


rimax = length(result(1).rr)

dy =5;


imgiR = zeros(rimax*dy,timax);
imgiG = imgiR;

for ti = 1:timax
   for i = 1:dy-1
      imgiR(i:dy:(dy*(rimax-1)+i), ti) = result(ti).cRr;
      imgiG(i:dy:(dy*(rimax-1)+i), ti) = result(ti).cGr;
   end
end


% scale
imgiR = imgiR / (1.2 * 10^-3);
imgiG = imgiG / (1.2 * 10^-3);

figure(111); clf;
%imsubplot(2,1,1);
%implot(imgray2color(imgiG', 'g'))
%imsubplot(2,1,2);
%implot(imgray2color(imgiR', 'r'))


imgi = cat(2, imgray2color(imgiG', 'g'), imgray2color(imgiR', 'r'));
implot(imgi)


%%
imwrite(impqlpermute(imgi, 'pqc', 'qpc'), fullfile(exp.ResultDirectory, ['rosette_class_evolution.tif']))



%% Intensities as Density Plots

ti = 1;

for ti = 1:timax

   ti
   
   isize = size(imglab);
   rmax = 700; rmin = 100;
   dr = 50;
   rr = rmin:dr:rmax;
   rmax = rmax + dr;
   
   th1 = 3/4*pi - 0.1;
   th2 = 5/4 * pi;

   ri = 1;
   clear cR cG
   cR= zeros(1, length(rr));
   cG = cR;
     
   stats = fstats{ti};
   
   for r = rr
      %count green / red
      
      id = and([stats.dist] >= r, [stats.dist] < r+dr);
      id = and(id, and([stats.theta] >= th1 , [stats.theta] <= th2));

      idR = and(id, [stats.class] == 1);
      idG = and(id, [stats.class] == 2); 
      
      cR(ri) = median([stats(id).RIntensity]);
      cG(ri) = median([stats(id).GIntensity]);

      %{sum(id), cR(ri), cG(ri)}
      
      %ar(ri) = overlapDiskSegmentRectangle(rcenter, r + dr, th1, th2, [0,0], isize) - overlapDiskSegmentRectangle(rcenter, r, th1, th2, [0,0], isize);
      
      %cRr(ri) = cR(ri) / ar(ri);
      %cGr(ri) = cG(ri) / ar(ri);
      
      ri = ri +1;
   end
   
   ires(ti).cR = cR;
   ires(ti).cG = cG;
   
end


%%

rimax = length(rr)

dy =5;

imgiR = zeros(rimax*dy,timax);
imgiG = imgiR;
size(imgi)

for ti = 1:timax
   for i = 1:dy-1
      imgiR(i:dy:(dy*(rimax-1)+i), ti) = ires(ti).cR;
      imgiG(i:dy:(dy*(rimax-1)+i), ti) = ires(ti).cG;
   end
end

{max(imgiR(:)), max(imgiG(:))}

%%
% scale
imgiR = imgiR / 0.55;
imgiG = imgiG / 0.55;

figure(111); clf;
%imsubplot(2,1,1);
%implot(imgray2color(imgiG', 'g'))
%imsubplot(2,1,2);
%implot(imgray2color(imgiR', 'r'))


imgi = cat(2, imgray2color(imgiG', 'g'), imgray2color(imgiR', 'r'));
implot(imgi)

%%

imwrite(impqlpermute(imgi, 'pqc', 'qpc'), fullfile(exp.ResultDirectory, ['rosette_intensity_evolution.tif']))






%% Plot Profile Time Evolution

rr  = result(1).rr;

for i = 1:length(result)
   figure(54); clf;

   %subplot(2,1,1)
   hold on
   plot(rr, result(i).cRr, 'r*-')
   plot(rr, result(i).cGr, 'g*-')

   ylim([0, 120*10^-5])
   xlim([0, rmax])
   
   pause(0.4)

   %subplot(2,1,2)
   %hold on
   %plot(rr, cR, 'r*-')
   %plot(rr, cG, 'g*-')

   %subplot(3,1,3)
   %plot(rr, ar, 'b')
end


%% Plot Profiles

rr  = 100:100:700;
dr = 100;
th1 = 3/4*pi - 0.1;
th2 = 5/4 * pi;

rc = exp.loadData('rosette_centers.mat');
isize = [ 1344,    1024 ];

tmax = length(result);
t2 = 100;

for ti = [1,t2, tmax]
   
   stats = fstats{ti};
   
   ri = 1;
   clear cR cG cRr cGr ar
   cR= zeros(1, length(rr));
   cG = cR; cRr = cR; cGr = cR; ar = cR;
   
   for r = rr
      %count green / red
      
      id = and([stats.dist] >= r, [stats.dist] < r+dr);
      id = and(id, and([stats.theta] >= th1 , [stats.theta] <= th2));
           
      cR(ri) = sum([stats(id).class] == 1);
      cG(ri) = sum([stats(id).class] == 2);
      
      %{sum(id), cR(ri), cG(ri)}
      
      % Rosette Center
      rcenter = fix(rc(ti, (1:2)));
      
      ar(ri) = overlapDiskSegmentRectangle(rcenter, r + dr, th1, th2, [0,0], isize) - overlapDiskSegmentRectangle(rcenter, r, th1, th2, [0,0], isize);
      
      cRr(ri) = cR(ri) / ar(ri);
      cGr(ri) = cG(ri) / ar(ri);
      
      ri = ri +1;
   end
   

   h = figure(200+ti); clf;

   %subplot(2,1,1)
   hold on
   plot(rr, cRr, 'r*-')
   plot(rr, cGr, 'g*-')

   ylim([0, 120*10^-5])
   xlim([0, rmax])
   
   pause(0.4)

   %subplot(2,1,2)
   %hold on
   %plot(rr, cR, 'r*-')
   %plot(rr, cG, 'g*-')

   %subplot(3,1,3)
   %plot(rr, ar, 'b')
   
   
   print(h,'-dpdf', exp.ResultFile(['resotte_radial_density_T' num2str0(ti, 3) '.pdf']))
   
   
end





%% Remove low Intensity Cells


figure(9);clf
subplot(1,2,1)
hist(rdat)
subplot(1,2,2)
hist(gdat)


% blue cells

id = rdat + gdat < 0.325;
%id = and(id, gdat < 0.25);
id = or(id, ddat > 850);
cdat3 = cdat;
cdat3(id) = 3;


%
h = figure(17); clf;
hold on

id = (cdat3 == 2);
plot(tdat(id), ddat(id), 'g.');
id = (cdat3 == 1);
plot(tdat(id), ddat(id), 'r.');

xlim([0,215]); ylim([0, 900])


%%
exp.saveData('mathematica.mat', [tdat; ddat; rdat; gdat; cdat; cdat3], '-v4')

%%
print(h,'-dpdf', exp.ResultFile(['resotte_radius_cells.pdf']))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Scaling



scalemin = [1,1,1]*2^16;
scalemax = [1,1,1]*0;

tstart = 7;
tend = 217;
%tend = 7;
tdelta = 8;
timax = length(tstart:tdelta:tend);


for t = times(tstart:tdelta:tend)

   fprintf('time: %g  index: %d / %d\n', t, ti, timax); 
   
   % load image
   imgB = imread(exp.fileName('time', t, 'channel', channelnames{1}))';
   imgR = imread(exp.fileName('time', t, 'channel', channelnames{2}))';
   imgG = imread(exp.fileName('time', t, 'channel', channelnames{3}))';

   img(:,:,1) = imgR; img(:,:,2) = imgG; img(:,:,3) = imgB;

   for c = 1:3
      scalemin(c) = min(scalemin(c), min(min(img(:,:,c))));
      scalemax(c) = max(scalemax(c), max(max(img(:,:,c))));
   end
end

scalemin
scalemax



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Rosette Center


tt = 1:15:217

timax = length(tt);

imgall = zeros(1344, 1024, timax, 3);

ti = 1;
for t = tt
   imgall(:,:,ti,:) = fucci_load(exp, t);
   ti = ti+1;
end


cent = zeros(timax, 2);

for ti = 1:timax
   figure(53); clf;
   implot(squeeze(imgall(:,:,ti,:)))
   drawnow
   
   pt = impoint;
   
   fprintf('center for time %d is %s\n', tt(ti), var2char(pt.getPosition()));
   
   cent(ti, :) = fix(pt.getPosition());
   
   %pause(0.4)
end


%%

for ti = 1:timax
   centt(ti,:) = [cent(ti,:), tt(ti)];
end
centt

exp.saveData('rosette_centers_raw.mat', centt)

%% interpolated center

centt = exp.loadData('rosette_centers_raw.mat')

t = 1:217;

cx = interp1(centt(:,3), centt(:,1), t, 'linear','extrap');
cy = interp1(centt(:,3), centt(:,2), t, 'linear','extrap');

figure(51); clf;
plot(cx,cy)

centall = [cx',cy']

exp.saveData('rosette_centers.mat', centall)



%% Check Centers

rc = exp.loadData('rosette_centers.mat');
t = 1;

%%
t
img = fucci_load(exp, t);

figure(1); clf;
implot(img); hold on
plot(rc(t,1), rc(t,2), '*')

t = t+1;


%% Automatic Rosette Center Finder - old

%rcenter = fucci_center(img, verbose);
%rcenter = rcenter'






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Segmentation and Classification

% use nuclear + green + red
t = 201;

img = fucci_load(exp, t);

img = 1 * img(:,:,1) + 1*img(:,:,2) + 1 * img(:,:,3);

figure(1); clf; imcolormap('b')
implot(img)


%%
%imgmask = img > 0.125;
imgmask = img > 0.25;
imgmask = imopen(imgmask, strel('disk', 3));
imgmask = postProcessSegments(bwlabeln(imgmask), 'volume.min', 50) > 0;

if verbose
   max(img(:))
   
   figure(21); clf;
   set(gcf, 'Name', ['Masking'])
   implot(imoverlaylabel(img, double(imgmask)))
end

%% Seeds

imgf =img;
imgf = filterLoG(max(imgf(:)) - imgf, [32,32], []);
%imgf = mat2gray(imgf);

{min(imgf(:)), max(imgf(:))}

% h-max detection (only local maxima with height > hmax are considered as maxima
imgmax = imextendedmax(imgf,  0.00005);

%local max
%imgmax = imregionalmax(imgf);

% constrain to maxima within mask
imgmax = immask(imgmax, imgmask);

% Combine nearby points
imgmax = imdilate(imgmax, strel('disk', 4));
% fill holes - combination of nearby points can lead to holes
imgmax = imfill(imgmax,'holes');
% shrink to single points - extended maxima usually give better segmentation results
% imgmax = bwmorph(imgmax,'shrink',inf);           

% plot the results.
if verbose  
   figure(22); clf
   set(gcf, 'Name', ['Seeding'])
   implottiling({imoverlay(img, imgmax), imoverlay(imgf, imgmax)});
end


























%% Plot as Dots


figure(50); clf; 

%imsubplot(1,2,1)
%imcolormap(gray)
imgp = squeeze(imgall(:,:,ti,:));
implot(imgp);  hold on
plot([rcenter(2)], [rcenter(1)],'Color','r','LineWidth',2, 'Marker', '*', 'MarkerSize', 20)
%rectangle('Position',[rcenter,50,50], 'Curvature',[1,1], 'FaceColor','r')

%subplot(1,2,2)
hold on
posR = pos(:, class==1);
posG = pos(:, class==2);

plot(posR(1,:),posR(2,:), '.r')
plot(posG(1,:),posG(2,:), '.g')

plot([rcenter(1)], [rcenter(2)],'Color','b','LineWidth',2, 'Marker', '*', 'MarkerSize', 20)




