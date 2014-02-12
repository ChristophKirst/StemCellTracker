%% fake image


image = fspecial('disk', 6);
image = mat2gray(padarray(image, [10 10]));

imagegrad = imgradient(image);


figure(1)
clf
ax(1) = imsubplot(1,2,1);
imshow(image);
ax(2) = imsubplot(1,2,2);
imshow(imagegrad);

linkaxes(ax, 'xy')


segmentByActiveRays(image, imagegrad, [], [] ,[])