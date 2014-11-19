
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Serial Detection of Colonies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data is read from disk each time and not cached
% useful for large tilings that do not fit into memory

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;

%initializeParallelProcessing(4)


tic

%% Overview Image

%isOverview = ImageSourceBF('/home/ckirst/Data/Science/Projects/StemCells/Experiment/Cytoo_IF/131106_RUES2_Lam521_46hBMP4_Bra_Sox2_Cdx2/Cy_w1_s3009_t1.TIF')

%figure(1)
%isOverview.plot

%% Setup Image Source
clc

% infer the tag expression o he files from the file folder automatically
%direc = '/Volumes/Gold digga/RUES2 experiments/11082014/LY40_DAPI_OCT4_NANOG_SOX2_01';
%texp = '/Volumes/Gold digga/RUES2 experiments/11082014/LY40_DAPI_OCT4_NANOG_SOX2_01/image_w435_s<S>_t1.TIF';


texp = '/home/ckirst/Science/Projects/StemCells/Data/TrackingPC/Anna/11082014/LY40_DAPI_OCT4_NANOG_SOX2_01/image_w435_s<S>_t1.TIF';

% fns = tagExpressionToFiles(texp);
% length(fns)
%load([direc filesep 'outall'],'dims')

%%


%%
clc
is = ImageSourceFiles(texp);
%is.fromFileExpression(texp, 'S', 1:2162);

% set the tiling and cell formats
%is.setReshape('S', 'UV', [46, 47]);
%is.setCellFormat('Vu');
is.printInfo


%%

% cd = is.cell('U', 1:4, 'V', 1:5);
% figure(2); clf
% implottiling(cd)


%% full preview
% 
is.addRange('U', 1:15, 'V', 1:15);
% 

preview = is.previewStiched('overlap', 120, 'scale', 0.05, 'lines', false);
% 
figure(2); clf
implot(preview)

%% restric range to some sub set

%is.addRange('U', 1:4, 'V', 1:4)
%is.setRawCellDataCaching(true);
%figure(1); clf
%is.plottiling
