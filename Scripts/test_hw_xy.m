image = fspecial3('disk', [20,20,20],5,1,0)
image = padarray(image, [20,10,0],'pre');

[x,y,z] = imfind(image)

size(image)


figure(1)
clf; imshow3d(image)

figure(2)
clf; implot3d(image)


