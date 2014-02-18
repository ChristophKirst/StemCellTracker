%% Test 3D segmentation of Arye Data


%% clear 

close all
clear all
clc

fp = [1   705   560   420];
set(0, 'DefaultFigurePosition', fp)

initialize()

%% other test data - "intensity"

x0 = 200; y0 = 150; wh = 150;
stack = imreadstack('./Test/Images/Stack2/*c001.tif', 'PixelRegion', {[x0 x0+wh], [y0 y0+wh]});
stack = imzinvert(stack);
stack = imisostack(stack, 2);
size(stack)


boxr = [1, 1, 1];

figure(1)
clf
imshow3d(stack, 'BoxRatios', boxr )
xlabel('x'); ylabel('y'); zlabel('z');

%figure(2)
%immontage(stack)

figure(3)
hist(log2(stack(stack > 0)), 120)


%%ijplot3d(stack)


%% get everything from Imaris

stack = imarisget();
stack = imisostack(stack, 2);
boxr = [1, 1, 1];
size(stack)

figure(1)
clf
imshow3d(stack, 'BoxRatios', boxr )
xlabel('x'); ylabel('y'); zlabel('z');

%figure(2)
%immontage(stack)

%figure(3)
%hist(log2(stack(stack > 0)), 120)

%% threshold

th = 2^thresholdMixtureOfGaussians(log2(stack(stack>0)), 0.1);
stackth = stack;
stackth(stack < th) = 0;
stackmask = stackth > 0;
stackmed = medianFilter(stackth, [3 3 3]);


figure(10)
clf
imshow3d(stackth, 'BoxRatios', boxr )

figure(11)
clf
imshow3d(stackmed, 'BoxRatios', boxr )

%%

figure(12)
immontage(stackmed)


%% log and maxima

stacklog = logFilter(max(stackmed(:)) - stackmed, [7 7 7]);
%stacklog = imhmax(mat2gray(stacklog),0.1); 
stackmax = imextendedmax(mat2gray(stacklog), 0.05);

% 3d water shed 
rm = imimposemin(max(stackmed(:)) - stackmed, stackmax);
ws = watershed(rm);


figure(20)
clf
colormap(autumn)
is = stacklog- mean(stacklog(:));
is(is < 0) = 0;
imshow3d(is, 'BoxRatios', boxr);

figure(21)
clf
imshow3d(stackmax , 'BoxRatios', boxr )

figure(22)
clf
imshow3d(double(ws) .* double(stackmask))


%% postprocess

seg = double(ws) .* double(stackmask);

seg2 = bwlabeln(seg > 0);

length(unique(seg(:)))
length(unique(seg2(:)))


% segments size

area = regionprops(seg2, 'area');
area = [area.Area];

%px = regionprops(seg2, 'PixelIdxList');
%px = {px.PixelIdxList};

idx = find(area < 50);
segth = seg2;
for i = idx
   segth(segth == i) = 0;
end

% fill possible holes 

segh = segth;
for s = 1:size(segth,3)
        segh(:,:,s) = imfill(segh(:,:,s), 'holes');
end

% remove boundary objects

% we dont want to clear objects that border in z
pseg = zeros(size(segh) + [0 0 2]);
pseg(:, :, 2:end-1) = segh;
pseg = imclearborder(pseg);
segh = pseg(:,:,2:end-1);
segh = bwlabeln(segh>0);

props = regionprops(segh, stack, 'Area', 'MaxIntensity', 'MinIntensity','MeanIntensity', 'PixelValues');


% plotting 
figure(30)
clf
imshow3d(seg2);

figure(31)
hist(area, 50)

figure(32)
clf
imshow3d(segh)


figure(33)
nbin = 25;
subplot(2,3,1)
hist([props.Area], nbin)
title('Area [Pixel]')
subplot(2,3,2)
hist([props.MaxIntensity], nbin)
title('Max Intensity')
subplot(2,3,3)
hist([props.MinIntensity], nbin)
title('Min Intensity')
subplot(2,3,4)
hist([props.MeanIntensity], nbin)
title('Mean Intensity')


%% select some shapes

idx = unique(segh(:));
if ~isempty(idx) && idx(1) == 0
   idx = idx(2:end)';
end
idx

% remoce some points
idx = idx(50:55);

seghs = segh;
seghs(seghs < 50) = 0;
seghs(seghs > 55) = 0;
seghs = bwlabeln(seghs);


figure(40)
clf
implot3d(seghs);


%% plot as surfaces

clear props
i = 1;
for l = idx
   ind = find(segh == l);
   [x,y,z] = ind2sub(size(segh), ind);
   props(i).ConvexHull = convhulln([x,y,z]);
   props(i).PixelPositions = [x,y,z];
   i = i + 1;
end

figure(41)
clf
hold on
for i = 1:length(props)
   trisurf(props(i).ConvexHull, props(i).PixelPositions(:,1),props(i).PixelPositions(:,2), props(i).PixelPositions(:,3));
end


%% surface pixel

idx = unique(seghs(:))';
surf= cell(length(idx),1);
i = 1;
for l = idx;
   surf{i} = bwsurface(seghs == l);  
   i = i + 1;
end


figure(50)
clf
implot3d(surf{2});



%hist([props.MeanIntensity], nbin)
%title('Surface Area')


%% imaris

obj = seghs == 2 ;

imarisput(255 * uint8(imisostack(obj,-2)))


%%
clf
implot3d(obj)

%% set surface in imaris

[xx,ff,nn] = imisosurface(obj);

imarisput(xx,ff,nn);



%% what about sherical filters


%stacksp = sphereFilter(stackmed, [8, 8, 4]);
stacksp = diskFilter(stackmed, [8, 8, 6], 1, 1, -0.2);
stackmax = imextendedmax(mat2gray(stacksp), 0.0);

% 3d water shed 
rm = imimposemin(max(stackth(:)) - stackth, stackmax);
ws = watershed(rm);


figure(40)
clf
colormap(autumn)
is = stacksp- mean(stacksp(:));
is(is < 0) = 0;
imshow3d(is, 'BoxRatios', boxr);


figure(41)
clf
imshow3d(stackmax , 'BoxRatios', boxr)

figure(42)
clf
imshow3d(double(ws) .* double(stackmask))


figure(43)
clf
imshow3d(imoverlay

