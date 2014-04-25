%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageFormats Interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
initialize

%%
bfinitialize

%%%%%%%%%%%%%%%%%
%% imread_bf_info


%% tif 
meta = imread_bf_info('./Test/Images/hESCells_DAPI.tif')
meta.parameters


%% lif
meta = imread_bf_info('/home/ckirst/Science/Projects/StemCells/Experiment/Other/Ian/Zstack_40x_50um_stepsize_2um.lif')
meta.parameters;


%% lsm

lsmfile = '/home/ckirst/Media/ChristophsData/Science/Projects/StemCells/Experiment/Other/Aryeh/YFPCitrine_movie1.lsm'
meta = imread_bf_info(lsmfile)
nseries = length(meta)

%% 
img = imread_bf(lsmfile,struct('series', 1, 'time', 25));
size(img)
max(img(:))
min(img(:))
class(img)

%%
figure(1); clf;
implottiling(img/4095)

%%
for z=1:size(img,3)
   imwrite(img(:,:,z)/4095, ['./Test/Images/hESCells_tif_Citrine/img_z' num2str0(z, 3) '.jpg']);
end

