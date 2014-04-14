%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup an Experiment, Segment and Analyse %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

initialize

%% Generate Experiment Class 


exp = Experiment('Name', 'Example', 'Date', datestr(now), 'Description', 'Test the Code',...
                 'BaseDirectory', './Test/Data/Experiment', ...
                 'ImageDirectoryName', '../../Images/hESCells_tif_stack',...
                 'ResultDirectoryName', '.', ...
                 'ReadImageCommandFormat', 'imread(''<directory>/<file>'')',...
                 'FileFormat', 'W1F127T<time,4>Z<z,2>C1.tif', ...
                 'FileFormatNames', {'z', 'time'});

exp.Info()

%% Save/Load some test data

sdata = rand(1,3);

exp.SaveResult('testdata.mat', sdata)

ldata = exp.LoadResult('testdata.mat');
ldata - sdata

%% Load Images

exp.ReadImageCommand([4, 1])

img = exp.ReadImage([4, 1]);
figure(1); clf;
implot(img)




