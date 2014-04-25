figure(1); clf

img = uint8(255 * rand(10,20));

colormap jet

cm = colormap;
imshow(img, 'DisplayRange', [0,max(img(:))])
%colormap(cm)


%%
figure(2); clf
img = uint8(255 * rand(10,20));

img = rand(30,24);

colormap jet
imsubplot(1,3,1);
implot(img)

colormap gray
imsubplot(1,3,2);
implot(img, 'color.scale', [0, 0.2])

imgc = imgray2color(img, 'red');
imsubplot(1,3,3)
implot(imgc)

colormap jet