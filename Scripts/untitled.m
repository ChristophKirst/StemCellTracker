%% Testing 3D Segmentation

path(path, './Test');


%% synthetic data
[img, ctr] = syntheticConfocalData();

figure(42)
vol3d('cdata', img)

figure(43)


img3D = zeros(size(img,1), size(img,2), 1, size(img,3));
for z=1:size(img,3)
   img3D(:,:,1,z) = mat2gray(img(:,:,z));
end

montage(img3D)


%% test data


dir