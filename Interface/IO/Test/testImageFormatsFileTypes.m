%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageFormats File types %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
initialize
bfinitialize


%% tif 
ii = imreadBFInfo('./Test/Images/hESCells_DAPI.tif')

ii.imetadata


%% tif stack 

ii = imreadBFInfo('./Test/Images/hESCells_tripple_reporter.TIF')


%% lsm
 
ii = imreadBFInfo('./Test/Images/hESCells_100x.lsm')

ii.imetadata

%% larger lsm

ii = imreadBFInfo('./Test/Images/hESCells_DAPI_Phallodin_ECad_CDX2_63x.lsm')


%% avi

ii = imreadBFInfo('./Test/Images/hESCells_tripple_stack.avi')

ii.imetadata


%% lif file




