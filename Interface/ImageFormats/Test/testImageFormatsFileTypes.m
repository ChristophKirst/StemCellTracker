%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageFormats File types %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all
initialize
bfinitialize


%% tif 
ii = imread_bf_info('./Test/Images/hESCells_DAPI.tif')

ii.imetadata


%% tif stack 

ii = imread_bf_info('./Test/Images/hESCells_tripple_reporter.TIF')


%% lsm
 
ii = imread_bf_info('./Test/Images/hESCells_100x.lsm')

ii.imetadata

%% larger lsm

ii = imread_bf_info('./Test/Images/hESCells_DAPI_Phallodin_ECad_CDX2_63x.lsm')


%% avi

ii = imread_bf_info('./Test/Images/hESCells_tripple_stack.avi')

ii.imetadata


%% lif file




