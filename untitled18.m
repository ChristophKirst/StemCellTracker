figure(1); clf

img = uint8(255 * rand(10,20));

colormap jet

cm = colormap;
imshow(img, 'DisplayRange', [0,max(img(:))])
colormap(cm)