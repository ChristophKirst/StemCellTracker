%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FileHandler and Experiment class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of the FileHandler class is to organize file access and writing
% so that reading image data eventually becomes indipendent of the underlying data format obtained from
% various microscopes etc...

clear all
clear classes
close all
clc

initialize


%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FileHandler Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

fh = FileHandler('BaseDirectory',          './Test/Data/Experiment', ...
                 'ImageDirectoryName',     '../../Images/hESCells_Cytoo',...
                 'ResultDirectoryName',    '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat',    'Big1_CFP_<time>.tif');

fh.info()

%%

fh.ReadImageTagNames


%%

fns = fh.files'


%%

fh.ReadImageCommand('time', 1)
fh.fileName('time', 1)
fh.fileName

%%

img = fh.readImage('time', 3);
figure(42); clf; colormap gray
implot(img)


%%
img = fh.readImage('time', 1);
img = mat2gray(img);

figure(1)
colormap(gray)
implot(img)

figure(2)
subplot(1,2,1);
hist(img(:), 256)
subplot(1,2,2);
hist(log(img(:))+eps, 256)


%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% The Experiment class extends the FildHandler class by adding information usefull for the experiment,
% it also stores the results


clear all
close all
clear classes
clc

%% Generate Experiment

exp = Experiment('name', 'Example', 'date', datestr(now), 'description', 'Test the Code',...
                 'BaseDirectory', './Test/Data/Experiment', ...
                 'ImageDirectoryName', '../../Images/hESCells_Cytoo',...
                 'ResultDirectoryName', '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat', 'Big1_CFP_<time>.tif');

exp.info()

%% Save/Load some test data
 
sdata = rand(1,3);

exp.saveData('testdata.mat', sdata)

ldata = exp.loadData('testdata.mat');
ldata - sdata

%% Load Images and plot

img = exp.readImage('time', 2);

figure(1); clf;
implot(img)


%% Save the Experiment 

imglab = syntheticLabeledImage([40, 40], 5);
img = rand(40,40);

objs = label2DataObjects(imglab, img);

exp.result = objs;


exp.saveExperiment('experiment.mat')

%%
clear all

%%
exp = loadExperiment('./Test/Data/Experiment/experiment.mat')

objs = exp.result;
[objs.r]
[objs.id]


















