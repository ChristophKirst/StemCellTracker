%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageInfo class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
clear classes
close all

initialize
bfinitialize


%% empty class

ii = ImageInfo

ii.size
ii.dim

%% from image file
 
ii = imread_bf_info('./Test/Images/hESCells_DAPI.tif')

ii.size
ii.dim                                              

%%

ii.setSize([100,200])


%% Test reformatting

ii.renameInFormat('p','c')

%%

ii.permuteToFormat('qc')

%%

ii.formatpos('c')
ii.dataformatpos('c')

%% series integration

ii = ImageInfo


ii.isize = [100,200,3];
ii.iformat = 'pqc';

ii.iseries = 1:100;
ii.iseriesdim = 'c';


ii.pqlctsizeFromFormatAndSize


ii.seriesAssignmentIds



