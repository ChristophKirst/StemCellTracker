%% Test Triangulation Tools


%% Create Test Object

image = fspecial3('disk', [40,20,15], [1, 2, 2], 1, 0);
image = padarray(image, [20,10,0], 'pre');

[x,y,z] = imfind(image);
size(image)


%figure(1)
%clf; imshow3d(image)

figure(2)
clf; implot3d(image)

line([1,40], [1,30],[1,20])

obj = image;

%% imaris

iim = imarisget();
size(iim)

obj = zeros(size(iim));
obj = imreplace(obj, image, [400, 70, 2]);

figure(1)
clf; implot3d(obj)

figure(2)
clf; implot3d(iim)


%% check if we get smae 3d picture in imaris

imarisput(205 * uint8(obj))
% ok, works

%% imisosurface

[x,y,z] = imfind(obj);
tri = delaunayn([x,y,z]);

figure(50)
trisurf(tri,x,y,z);
xlabel('h'); ylabel('w');

[xyz, tri, nrm] = imsurface(obj);

figure(3)
clf; imsurfaceplot3d(xyz, tri, nrm)



%% put surface into imaris 

imarisput(xyz, tri, nrm)



%% 

istack = imarisget();
size(istack)

label = istack > 0;

[v, f, nrm] = imisosurface(label);


imarisput

