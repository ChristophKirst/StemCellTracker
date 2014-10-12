%%%%%%%%%%%%%%%%%%%%
%%%  Test Filter %%%
%%%%%%%%%%%%%%%%%%%%

initialize

%% 2D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Test Image

img = imread('./Test/Images/hESCells_DAPI.tif');
img = mat2gray(img);

figure(1)
implot(img)

%% Gaussian Filter

gf = fspecial2('Gaussian', [15, 8],[3,8]);
total(gf)

figure(2)
imagesc(gf)
colormap jet
colorbar

imgf = filterGaussian(img, 15, 15);

figure(3)
implottiling({img, imgf})


%% LoG Filter

gf = fspecial2('LoG', [15, 15],[4,4]);
total(gf)
{max(gf(:)), min(gf(:))}

figure(2)
imagesc(gf)
colormap jet
colorbar


imgf = filterGaussian(img, 15, 15);

figure(3)
implottiling({img, imgf})


%% %% DoG Filter

gf = fspecial2('DoG', [15, 15],[4,4]);
total(gf)
{max(gf(:)), min(gf(:))}

figure(2)
imagesc(gf)
colormap jet
colorbar


imgf = filterGaussian(img, 15, 15);

figure(3)
implottiling({img, imgf})



%% Sphere Filter

gf = fspecial2('Sphere', [15, 15],[7,17]);
total(gf)
{max(gf(:)), min(gf(:))}

figure(2)
imagesc(gf)
colormap jet
colorbar


imgf = filterGaussian(img, 15, 15);

figure(3)
implottiling({img, imgf})



%% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


