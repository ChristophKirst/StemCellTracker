%% test the Bio Fromats interface

clear all
close all
clc

initialize
bfinitialize


%% load image via file browser

img = imread_bf()

%% load image data from a lsm file

data = imread_bf('/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Other/Aryeh/YFPCitrine_movie1.lsm', struct('series', 1, 'time', 1));
size(data)


%% load image data from a lif file

