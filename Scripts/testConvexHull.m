close all
clear all


%% test 3d convex hull

mask = zeros(5,5,5);
mask(2,3,2) = 1;
mask(2,3,3) = 1;
mask(2,5,3) = 1;
mask(5,5,2) = 1;
mask(2,2,2) = 1; 
mask(5,2,5) = 1;



figure(1)
clf
imshow3d(mask)

siz = size(mask);

ind = find(mask);
[x,y,z] = ind2sub(siz, ind);


[gx, gy, gz] = meshgrid(1:siz(1), 1:siz(2), 1:siz(3));

[A,b] = vert2con([x, y, z]);

r = [gx(:)' ; gy(:)' ; gz(:)'];

id = all(A * r - repmat(b,1, length(r))<= 0.5);

conv = zeros(siz);
conv(id > 0) = 1;

figure(2)
clf
imshow3d(conv);