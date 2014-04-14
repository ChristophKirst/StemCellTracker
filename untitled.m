%% test labeling etc


imglab = syntheticLabeledImage([50, 50], 5);

figure(1); clf
implot(imcolorize(imglab))


%%

cent  = imstatistics(imglab, 'Centroid');
cent = round([cent.Centroid]);


isize = size(imglab);

clabel = zeros(isize);
clabel(imsub2ind(isize, cent')) = 1;


figure(2)
implot(imoverlay(imglab, clabel))


figure(3)
cl = bwlabeln(clabel);
implot(imcolorize(cl))