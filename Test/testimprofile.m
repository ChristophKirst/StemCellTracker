% Test coordinates for improfile

x = 1:20;
y = 1:30;
[x, y] = meshgrid(x,y);

z = (x-10).^2.* y;
zt = z';
zt = imreplace(zt, -2000 * ones(4,5), [13, 20]);
z = zt';


pos = find(zt == -2000);
[ppx, ppy] = ind2sub(size(zt), pos)


size(z)

px = [4, 19];
py = [10,25];


figure(1)
imagesc(mat2gray(z))
set(gca, 'YDir', 'normal')
line(px, py)


figure(2)
imshow(mat2gray(z))
set(gca, 'YDir', 'normal')
line(px, py)


figure(3)
implot(mat2gray(zt))
set(gca, 'YDir', 'normal')
line(px, py)

line(ppx, ppy)

z(py(1), px(1))
zt(px(1), py(1))


prof  = improfile(z, px, py, 20);

proft = improfile(zt, py, px, 20);
prof2 = improfile(zt, ppy, ppx);


figure(5)
clf
hold on
plot(prof)
plot(proft, 'r')
plot(prof2, 'g')


% -> conclusion: transposed coordinates are the pixel coordinates

