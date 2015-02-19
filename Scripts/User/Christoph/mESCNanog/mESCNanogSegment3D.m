function imgseg = mESCNanogSegment3D(imgraw, imgmask, verbose, time)

%% settings

tiling = [5,8];


%% Filter

 imgmed = filterMedian(imgraw, [5, 5, 3]);
 
 %imgmed = immask(imgmed, imgmask);
   
 if verbose
    figure(4)
    set(gcf, 'Name', ['Time: ', num2str(time), ' Masked Median']);
    implottiling(imgmed, 'tiling', tiling);
 end
  
 %%
 
   
 %%
 
 imgmask = mat2gray(imgf) > 0.05;
 
 if verbose
    figure(15); clf
    set(gcf, 'Name', ['Time: ', num2str(time), ' Mask']);
    implottiling(imoverlay(mat2gray(imgf(:,:,10:19)), imgmask(:,:,10:19), 'r', false), 'tiling', tiling);
 end
  
 %%
 
 
 imgfm = immask(imgf, imgmask);
 
 %%
 imgfmm = filterMedian(imgfm, [5,5,3]);
 
 
 
 
 %%
 
imggrad = zeros(size(imgfm));
 for z = 1:size(imgfmm,3)
  imggrad(:,:,z) = imgradient(imgfmm(:,:,z));
 end
 
 imggrad = mat2gray(imclip(imggrad, 0, 25));
 
 
%%
if verbose
    figure(15); clf
    set(gcf, 'Name', ['Time: ', num2str(time), ' Mask']);
    implottiling(imggrad(:,:,10:19), 'tiling', tiling);
end
  
 %%
 
 imgfmg  = imgfm - 0.5 * imclip(imggrad)
 
 
 if verbose
    figure(15); clf
    set(gcf, 'Name', ['Time: ', num2str(time), ' Mask']);
    implottiling(imggrad(:,:,10:19), 'tiling', tiling);
end
 
   %%
   
%    imgmedc = mat2gray(imgmed);
%    imgmedc(imgmedc  < 0.4) = 0; 
%    
%    if verbose
%       figure(4)
%       set(gcf, 'Name', ['Time: ', num2str(time), ' Masked Median']);
%       implottiling(imgmedc, 'tiling', tiling);
%    end
   
   
%% Seeding
%
% % improve seeding by adding gradients
%
% imggrad = zeros(size(imgmed));
% for z = 1:size(imgmed,3)
% imggrad(:,:,z) = imgradient(imgmed(:,:,z));
% end
%
% %figure(103)
% %clf
% %implot3d(stackgrad)
%
% mi = max(imgmed(:));
% mg = max(imggrad(:));
% stackmix = imgmed - 0.0 * mi / mg * imggrad;
% stackmix(stackmix < 0) = 0;
%
% if verbose
%
% figure(5)
% clf
% set(gcf, 'Name', ['Mixture: ' filenames]);
% implot3d(stackmix);
% %ijplot3d(stackmix)
%
% end
%%
% stacklog = filterLoG(iminvert(stackmix), [10,10,10]);
%
% figure(6)
% set(gcf, 'Name', ['filterLoG: ' filenames])
% clf
% implot3d(imclip(stacklog - mean(stacklog(:))))
%%
% stackdog = filterDoG(stackmix, [10,10,10]);
%
% figure(7)
% set(gcf, 'Name', ['filterDoG: ' filenames])
% clf
% implot3d(imclip(stackdog - mean(stackdog(:))))
%%
% stackdisk = filterSphere(stackmix, 10);
%
% figure(8)
% set(gcf, 'Name', ['filterDisk: ' filenames])
% clf
% implot3d(imclip(stackdisk - mean(stackdisk(:))))


%% select fitlered stack and detect seeds
%stackf = imclip(stackdisk - mean(stackdisk(:)));
%stackf = stackf .* cast(stackmask, class(stackf));

imgf = imgmed;
imgf = filterSphere(imgf, [16, 16, 6]);
%imgf = log(mat2gray(imgf)+eps) + 5;
%imgf(imgf <0) = 0;
imgf = mat2gray(imgf);


%figure(4); clf;
%implottiling(imgf, 'tiling', tiling);

%%
imgmax = imextendedmax(mat2gray(imgf), 0.015);
%stackmax = imregionalmax(stackf);
%imgmax = immask(imgmax, imgmask);
%
%stackmax = imerode(stackmax, strel('disk',3));
%stackmax = imdilate(stackmax, strel('disk',3));

if verbose
   % figure(10)
   % clf
   % set(gcf, 'Name', ['PreSeeding: ' filenames]);
   % is = imgf - mean(imgf(:));
   % is(is < 0) = 0;
   % implot3d(is);
   %
   % figure(11)
   % clf
   % set(gcf, 'Name', ['Seeding: ' filenames]);
   % implot3d(mat2gray(imgmax));
   % overlay in imagej
   stackc = gray2rgb(imgmax);
   %size(stackc);
   stackc(:,:,:,2:3) = 0;
   stackc = stackc + gray2rgb(mat2gray(imgraw));
   stackc(stackc > 1) = 1;
   
   figure(15); clf
   set(gcf, 'Name', ['Time: ', num2str(time), ' Seeding']);
   implottiling(stackc, 'tiling', tiling);
end


%% Segment

% 3d water shed on image + gradient
%mi = max(imgmed(:));
%mg = max(imggrad(:));
%imgws = imgmed - 0.25 * mi / mg * imggrad;
%imgws(imgws < 0) = 0;
imgws = imgmed;
imgwsmin = imimposemin(iminvert(imgws), imgmax);
imgseg = watershed(imgwsmin);
imgseg = immask(imgseg, imgmask);

if verbose
   %figure(13)
   %set(gcf, 'Name', ['Segmentation: ' filenames])
   %clf
   %implot3d(ws)
   %figure(14)
   %clf
   %set(gcf, 'Name', ['Segmentation: ' filenames])
   %imsurfaceplot3d(ws)
   %ijplot3d(imcolorize(ws) + gray2rgb(mat2gray(stackraw)))
   figure(14); clf
   set(gcf, 'Name', ['Time: ', num2str(time), ' Watershed']);
   imgovl = imoverlaylabel(imgmed, impixelsurface(imgseg), false);
   implottiling(imgovl, 'tiling', tiling);
end



end
