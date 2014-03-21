bfinitialize()

filename = '\\tracking-pc.rockefeller.edu\DATA\Aryeh\hESC\Cytoo_IF\Expt24_Snail_Bra_Sox2_timecourse\140305_RUES2_36hBMP4_Bra_Snail_Sox2.lif';

xr = [1000, 1300];  % use [] for all
yr = [1300, 1600];

cr = 4;
se = 21;
ti = 1;
img = imread_bf(filename, struct('series', se, 'time', ti, 'channel', cr, 'x', xr, 'y', yr));
img = imzreverse(squeeze(img));
   
   %%
   imNew=zeros(size(img));
   for ii=1:size(img,3)
   im2d=double(img(:,:,ii));
   imMed=medianFilter(im2d,6);
   %imNoPix=imerode(im2d,strel('disk',5));
   %imshow(imNoPix,[]);
    
    %im2d=medianFilter(im2d,3);
%    
%    
    thresh = 50;
    im2d(im2d>thresh)=imMed(im2d>thresh);
    im2d=medianFilter(im2d,3);
    imshow(im2d,[]);
    imNew(:,:,ii)=im2d;
   end
  
   
%%
figure(3)
clf
   
im3= imNew(:,:,3:end);
implot3d(im3)



%%

imarisinitialize
imarisstart

%%
dset = imarisgetdataset

%%

imarissetvolume(dset, uint8(im3),3)


%%



%% get statistics

seg = imarisgetobject('Segmentation')

istat = seg.GetStatistics()


%%
load('Z:\aryeh_sox2_snail.mat')


%%

id1 = [stats{3}.MeanIntensity] >= 14.4 & [stats{3}.MaxIntensity] >= 65.;

nset = size(id1);
sfset = surf(id1);
fcset = fac(id1);
nmset = norm(id1);

imarissetsurface('Snail_good', sfset, fcset, nmset);


%%

id2 = [stats{4}.MeanIntensity] >= 23.5 & ~id1;

nset = size(id2);
sfset = surf(id2);
fcset = fac(id2);
nmset = norm(id2);

imarissetsurface('Sox2_good', sfset, fcset, nmset);





%%
ids= ~id1 & ~id2;

nset = size(ids);
sfset = surf(ids);
fcset = fac(ids);
nmset = norm(ids);

imarissetsurface('None', sfset, fcset, nmset);













