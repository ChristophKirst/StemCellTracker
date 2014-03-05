%% Geodesic SKIZ and Segmentation

clear all
close all
clc

%% create test image

h = 70;w = 60;
img = zeros(h,w);
img = imreplace(img, fspecial2('disk', [16,16])> 0, [14,6]);
img = imreplace(zeros(h,w), fspecial2('disk', [40,20])> 0, [12,19]) + img > 0;
img = imreplace(zeros(h,w), fspecial2('disk', [10,20])> 0, [42,7]) + img > 0;
%img = imreplace(zeros(h,w), fspecial2('disk', [10,10])> 0, [36,5]) + img > 0;

mask = img > 0;

%mask = bwulterode(img);
seed = logical(index2mask([h,w], [21, 43, 47; 10, 30, 13]));

figure(1)
imshow(imoverlay(img, seed))


%% gedoesic distance

bwg = (1- mat2gray(bwdistgeodesic(img, seed)));

bwd = (1- mat2gray(bwdist(mask))) .* (img);

figure(2)

imsubplot(1,2,1)
imshow(bwg)

imsubplot(1,2,2)
imshow(bwd)


%% geodesic skiz = watershed on geodesic distance

label = bwlabel(seed);
skiz = imskiz(img, label);

figure(3)
imsubplot(1,2,1)
imshow(imoverlay(imcolorize(skiz), seed))
imsubplot(1,2,2)
imshow(imoverlay(img, seed))



%% create erosions

clear imge
imge{1} = img;
i = 1;
while any(imge{i}(:))
   imge{i+1} = imerode(imge{i}, strel('disk', 1));
   i = i + 1;
end
imge = imge(1:end-1);


figure(4)
imsubplot(1,3,1)

imgs = imge{1};
ne = length(imge);
for i = 2:ne
   imgs = imgs + imge{i};
end
imshow(imcolorize(imgs))
imsubplot(1,3,2)
ue = bwulterode(mask);
imshow(ue)

imsubplot(1,3,3)
imshow(imcolorize(skiz))


%% geodesic skiz of k th eroded image in k-1 eroded image

clear skiz
ne = length(imge);
skiz{1} = label;
for i = 2:ne
   skiz{i} = imskiz((imge{ne-i+1} + skiz{i-1}) > 0, skiz{i-1});
end


figure(13)
imsubplot(1,2,1)
imshow(imcolorize(skiz{end}))
imsubplot(1,2,2)
imshow(imcolorize(imgs))




