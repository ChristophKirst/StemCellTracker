%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wild type colony Detection %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clear class
close all
clc


initialize
bfinitialize

%%


texp = './Test/Images/hESCells_Tiling_WildType/140902_RUES2_BMP4_DAPI_GFP_R_Cy5__9_p<P,6>t00000001z001c<C,2,s>.tif';

is = ImageSourceFiles(texp);
is.printInfo



%%
clc
is.setReshape('P', 'UV', [8, 19]);
is.printInfo

is.cellSize


is.setRange('U', 1:3, 'V', 1:3, 'C', 1);
is.printInfo

is.cellSize


%%

is.setCellDataFormat('X', 'UvYC')

is.cellSize



%%
is.setCellDataFormat('XY', 'UvC');
is.printInfo

is.cellSize

%%
is.setRange('U', 1:3, 'V', 1:3, 'C', 1);
is.cellSize

%%

figure(1); clf
is.plottiling


%%

is.cell([1:6])