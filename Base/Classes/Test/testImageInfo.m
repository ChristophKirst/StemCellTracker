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

ii.dataSize
ii.dataSizeC


%% from image file
 
ii = imreadBFInfo('./Test/Images/hESCells_DAPI.tif')

ii.dataSize
ii.dataDims                                             


%% Test reformatting

ii.renameFormat('Y','Z')



%% test some reshape functionality
clc
ii.setReshape('X', 'XY', [16,32])



%% test with some data

d = rand(512, 512);

dr = ii.reshapeData(d);
size(dr)



%% ranges

clc

clear all
clear classes
close all


ii = imreadBFInfo('./Test/Images/hESCells_DAPI.tif')



%% 

clc
ii.setRange('X', 1:10);
ii.printInfo



