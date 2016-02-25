%%%%% Test 4D filter %%%%%

initialize


%%

basedir = '/home/ckirst/Desktop/test multi size segmentation/';
fns = {'560DAPI.png', '820DAPI.png', '1005DAPI.tif', '1324DAPI.png'};


%%


for i = 1:4
   
   imgraw{i} = double(imread(fullfile(basedir, fns{i})));
   
   figure(i); clf
   implot(imgraw{i})
   
   imgraw{i} = imgraw{i}(100:1000, 100:1000);

end

%%

for i = 1:4
   imgi{i}= iminvert(imgraw{i})
   
   figure(i); clf
   implot(imgi{i})
end






%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

for i = 1:4

%imgI = imgf;
%imgI = imgI - imopen(imgI, strel('disk', 20));
%imgI = mat2gray(imgI);


%imgf = filterMedian(imgraw, 2);
%imgf = imgraw;
%imgf = filterBM(imgraw, 'profile', 'np', 'sigma', 10);

imgf{i} = imgraw{i};

%imgg{i} = imgradient(imgf{i});
%imgf{i} = imgf{i} - 1 * imgg{i};
%imgff = imgg{i};
%imgff(imgff<0) = 0;
%imgf{i} = imgff;

%imgf = imgf - imopen(imgf, strel('disk', 20));

figure(1+i); clf
imsubplot(3,1,1)
implot(imgraw{i})
imsubplot(3,1,2)
implot(imgf{i})
imsubplot(3,1,3)
hist(imgf{i}(:), 256)

end



%%

for i = 1:4
   th{i} = thresholdFirstMin(imgf{i},'nbins', 50, 'delta', 20)
end


%%

imgWS = watershed(filterMedian(imgi{3}, 3));
figure(7); clf;
implottiling({imgraw{3}, imoverlaylabel(imgraw{3}, impixelsurface(imgWS))})



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:4

imgmask{i} = imgf{i} > th{i};
imgmask{i} = imopen(imgmask{i}, strel('disk', 4));


if verbose
   %max(img(:))
   
   figure(21+i); clf;
   colormap gray
   set(gcf, 'Name', ['Masking'])
   implottiling({mat2gray(imgf{i}); imgmask{i}}) 
end

imgfm{i} = immask(mat2gray(imgf{i}), imgmask{i});

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for d =3*( 10:5:50)
   imgfF = filterLoG(imgi{1}, d);

   imgfF = mat2gray(imgfF);

   imgfM = imextendedmax(imgfF, 0.05);
   
   imgfL = bwlabel(imgfM);
   
   nM(d) = max(imgfL(:));
   
   figure(60)
   implottiling({imgfF, imgfM, mat2gray(imgfL)}')
   drawnow
   
end


figure(33+i)
plot(nM)




%%

for i = 1:4
   imgg{i} = imgradient(imgf{i});
end

figure(7); clf
implottiling(imgg')


%%

for i = 1:4
   imgc{i} = xcorr2(imgg{i}, imgg{i});
end
figure(7); clf
implottiling(imgc')




%%
for i = [3,4]

imgM{i} = filterFunction(imgfm{i}, 8, 'min');

figure(7+i); clf;
implottiling({mat2gray(imgf{i}); imgM{i}})

figure(12+i); clf;
hist(imgM{i}(:), 256)

th = [0.05,0.05,0.05, 0.2];

imgmaskM{i} = imgM{i} > th(i);


figure(22+i); clf;
implot(imgmaskM{i})



nO = 20;
imgO = imgmaskM{i};
nLoss = zeros(1,nO+1);
nLoss(1) = sum(imgO(:));
for d = 1:50
   imgO = imopen(imgmaskM{i}, strel('disk', d));
   nLoss(d+1) = sum(imgO(:)); 
   
   %figure(60+d)
   %implot(imgO)
   
end
nLoss = diff(nLoss(end:-1:1));


figure(33+i)
plot(nLoss)

end









%%
for i = [2,3]

w = [10:5:60];
nw = length(w);

imgW{i} = zeros([size(imgf{i}), nw]);

for ww=1:nw
   imgW{i}(:,:,ww) = filterFunction(imgfm{i}, w(ww), 'min');
end

end


%%

for i = [2,3]

figure(7+i); clf;
implottiling(imgW{i}, 'tiling', [3,4])

end

%%





%% Seeds

imgBMI = sum(imgCf,3);
imgHSV = rgb2hsv(imgCf);

imgIf = imgHSV(:,:,3);
%imgIf = filterSphere(imgIf, 4);
%imgIf = filterDisk(imgIf, 8, 1, 1, -1);
imgIf = filterDoG(imgIf, 8);
imgIf = mat2gray(imgIf);
imgImax = imextendedmax(imgIf, 0.01);   
%imgImax = imregionalmax(imgIf);

stats = imstatistics(bwlabel(imgImax), {'PixelIdxList', 'MaxIntensity'}, imgIf);
statsF = stats([stats.MaxIntensity] < 0.35);
pxl = {statsF.PixelIdxList};
pxl = cat(1, pxl{:});
imgImax(pxl) = 0;

%figure(10); clf
%implottiling({imgC; imoverlay(imgIf, imgImax); imoverlay(imgC, imgImax, [1,1,1])})

imgmax = cell(1,nch);
for i = 1:nch
   imgmax{i} = imgCf(:,:,i);
   %imgmax{i} = filterSphere(imgmax{i}, 3); 
   %imgmax{i} = filterDisk(imgmax{i}, 8, 1, 1, -1);
   imgmax{i} = filterDoG(imgmax{i}, 8);
   imgmax{i} = mat2gray(imgmax{i});
   imgmax{i} = imextendedmax(imgmax{i}, 0.01);   
   
   stats = imstatistics(bwlabel(imgmax{i}), {'PixelIdxList', 'MinIntensity'}, imgCf(:,:,i));

   statsF = stats([stats.MinIntensity] < 0.1);
   pxl = {statsF.PixelIdxList};
   pxl = cat(1, pxl{:});
   imgmax{i}(pxl) = 0;
end

imgmaxAll = imgImax + imgmax{1} + imgmax{2} +0 * imgmax{3} > 0;

if verbose
   figure(6); clf
   implottiling({imoverlay(imgIf, imgImax), imoverlay(imgC, imgImax, [1,1,1]); ...
              imoverlay(imgCf(:,:,1), imgmax{1}),  imoverlay(imgCf(:,:,2), imgmax{2}); ...
              imoverlay(imgCf(:,:,3), imgmax{3}),  imoverlay(imgC, imgmaxAll, [1,1,1])}')
end

           
%% Watershed

imgWs = imimposemin(iminvert(imgIf), imgmaxAll);
imgWs = watershed(imgWs);
imgWs = immask(imgWs, imgmask);
imgWs = imlabelseparate(imgWs);


if verbose
   figure(7); clf;
   implot(imgWs)

   %imglab = bwlabel(imgWs);

   figure(8); clf;
   implottiling(imoverlay(imoverlaylabel(imgC, impixelsurface(imgWs), false), imgmaxAll, [1,1,1]))
end
   
imgS = imgWs;      


%% Postprocess

imgSP = postProcessSegments(imgS, imgIf, 'intensity.median.min', 0.025, 'volume.min', 15, 'fillholes', false);
imgSP = imrelabel(imgSP);

if verbose
   imgSp1 = impixelsurface(imgSP);
   imgSp1 = imoverlaylabel(imgCf, imgSp1, false);
   figure(5); clf;
   implot(imgSp1);
end

imglab = imgSP;

