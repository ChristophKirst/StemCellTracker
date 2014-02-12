%% test ImageJ

filename = './Test/Images/raw.tif';

% ImageJ watershed
ijread(filename);
MIJ.run('Threshold');

tic;
MIJ.run('Watershed');
toc

% Matlab watershed
ijread(filename);
MIJ.run('Threshold');
tsource = MIJ.getCurrentImage();
MIJ.run('Distance Map');
msource = MIJ.getCurrentImage();

tic;
ws = watershed(255-msource);
toc

dam = (ws==0)*255;
MIJ.createImage(or(dam, double(255-tsource))*255);


%% 3D plotting surfac

Z = membrane(1,25);

Imin = min(Z(:));
Imax = max(Z(:));
I = uint8( 200 * (Z-Imin) / (Imax-Imin) );

Miji(false);

imp = MIJ.createImage('MATLAB peaks', I, false);
%MIJ.createImage(image);

universe = ij3d.Image3DUniverse();
universe.show();

color = javax.vecmath.Color3f(240 / 255, 120 / 255, 20 / 255);
c = universe.addSurfacePlot(imp, ...
        javax.vecmath.Color3f(), ...
        'Matlab Peak in 3D', ...
        1, ...
        [true true true], ...
        1);

     
% universe.resetView();
% c.setColor(color);
% c.setTransform([1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]);
% universe.fireContentChanged(c);
% universe.centerSelected(c);
% universe.rotateUniverse(javax.vecmath.Vector3d(-1, -0.5, +0.2), +120 * pi / 180);   


%% 3D plotting volume


load('mri.mat');
I = squeeze(D);
[R, G, B] = ind2rgb(I, map);
R  = uint8(255 * R);
G  = uint8(255 * G);
B  = uint8(255 * B);
J = cat(4, R,G,B);
size(R)
size(J)


ijplot3d(J, 'PixelDepth', 2.5)


%% 3d plotting stack


stack  = imread_hdr('./Test/Images/Stack/StackIsoCropM');

figure(1)
clf
implot3d(stack)


ijplot3d(stack)




