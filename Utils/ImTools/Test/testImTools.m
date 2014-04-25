%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Tools %%%
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%
%% imforamt

img = rand(20,15,3);
imformat(img)


%%%%%%%%%%%%%%%%%%%%%%%%
%% implot2d

img = rand(15,20);
implot2d(img)

%implot uses pixel coordinates:
line([5 10], [10,11])


%%%%%%%%%%%%%%%%%%%%%%%%
%% implot3d

img = rand(15,20,4);

figure(1); clf;
implot3d(img)

%color
img = rand(15,20,4,3);
figure(2); clf;
implot3d(img)

% pixel coordinates
figure(3); clf;
img = zeros(15, 20, 5);
img(3,4,5) = 1;
img(10, 8, 2) = 1;

colormap jet
pp = implot3d(img)
line([3, 10], [4,8], [5,2])


%%%%%%%%%%%%%%%%%%%%%%%%
%% implot

img = rand(20,30);
figure(1); clf;
implot(img)

img = rand(20,30,3);
figure(2); clf;
implot(img)

img = rand(20,30,10);
figure(3); clf;
implot(img)


%%%%%%%%%%%%%%%%%%%%%%%%
%% imreplace

img = rand(20,30);
si = zeros(3,3);
imgr = imreplace(img, si, [4,4]);

figure(1)
implot(imgr)


%%%%%%%%%%%%%%%%%%%%%%%%
%% imcolor2rgb

rgb = imcolor2rgb('magenta')


%%%%%%%%%%%%%%%%%%%%%%%%
%% imgray2color

img = rand(20,30);
rgb = imgray2color(img, 'magenta');

figure(1)
implottiling({img, rgb})


%%%%%%%%%%%%%%%%%%%%%%%%
%% imcolorize

img  = syntheticLabeledImage([10, 20], 5, 4);
imgc = imcolorize(img);

colormap jet

cpre = colormap;
cpre(end, :)

figure(1); clf;
imshow(img);

cpost = colormap;
cpost(end, :)
%colormap is changed, but can even change it afterwards
colormap jet

%% 
cpre = colormap;
cpre(end, :)

figure(2); clf;
implot(img)

cpost = colormap;
cpost(end, :)

%colormap is fixed, active color map has not changed !
colormap jet

%% 
figure(1); clf
imsubplot(1,2,1)
implot(img)
imsubplot(1,2,2)
implot(imgc)
%different appropiate color maps per image


%%%%%%%%%%%%%%%%%%%%%%%%
%% imformat

img = rand(20, 30, 10, 3);
imformat(img)


%%%%%%%%%%%%%%%%%%%%%%%%
%% impqlpermute

img = rand(20, 30, 10, 3);
imgr = impqlpermute(img, 'pqlc', 'clpq');
size(imgr)

img(3,4,5,2)
imgr(2,5,3,4)

img = rand(20,30,1,3);
imgr = impqlpermute(img, 'pqlc', 'cqp');
size(imgr)

img(3,4,1,2)
imgr(2,4,3)

img = rand(20,30, 3);
imgr = impqlpermute(img, 'xyc', 'yxc');
size(imgr)
img(3,4,2)
imgr(4,3,2)

img = rand(2,4,3);
imgr = impqlpermute(img, 'xyc', 'pqc');
size(imgr)
img(2,4,2)
imgr(2,end-4+1,2)


%%%%%%%%%%%%%%%%%%%%%%%%
%% impqlpermute & implot, imshow

img = imcolorize(syntheticLabeledImage([20, 30], 5, 10));

figure(1); clf;
implot(img)



%%%%%%%%%%%%%%%%%%%%%%%%
%% imcolormap
cm = imcolormap(5)

cm = imcolormap('jet', 70);
size(cm)

img = rand(20,20);
%%
figure(1); clf
imsubplot(1,2,1)
cm = imcolormap('gray');
implot(img,[0,1])
imsubplot(1,2,2)
cm = imcolormap('igray')
implot(img,[0,1])



%%%%%%%%%%%%%%%%%%%%%%%%
%% imgray2color

img = rand(30,20);
imgc = imgray2color(img, 'g');
imgc2 = imgray2color(img, setParameter('color.map', 'jet'));

figure(1); clf
imsubplot(1,3,1)
colormap 'gray'
implot(img);
imsubplot(1,3,2)
implot(imgc);
imsubplot(1,3,3)
implot(imgc2);



%%

implottiling({img, imgc})




%%%%%%%%%%%%%%%%%%%%%%%%
%% imcolorize




%% Testing Others


figure(1); clf;
img = rand(20, 39);
imsubplot(1,2,1)
imshow(img*255, gray(255))
imsubplot(1,2,2)
imshow(img*255, jet(255))


figure(2); clf;
imsubplot(1,2,1)
implot(img, gray(255))
imsubplot(1,2,2)
implot(img, imcolormap('igray', 255))


figure(3); clf;
imsubplot(1,2,1)
imcolormap('gray')
implot(img)
imsubplot(1,2,2)
imcolormap('igray')
implot(img)
