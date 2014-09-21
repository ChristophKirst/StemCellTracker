clc
img = rand(1500,3000);

tic
rp = regionprops(logical(img), 'BoundingBox')
rp.BoundingBox
toc

tic
[a,b] = imboundingbox(img)
toc


