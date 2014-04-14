%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calcium Signalling in Stem Cells %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localroot = '/home/ckirst';
externalroot = '/var/run/media/ckirst/ChristophsData';
%rootdir = localroot;
rootdir = externalroot;

basedir = fullfile(externalroot, '/Science/Projects/StemCells/Experiment/Calcium');


%% GCamp6s Experiments

expdir = '/hES_GCamp6/2014_04_04';
dirname = fullfile(basedir, expdir);


info = imfinfo('/hES_GCamp6/2014_04_04/laser40_EM1000_exp100_001_t1.TIF')

%%
info = imfinfo('/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/laser40_EM1000_exp100_001_stack.TIF')

%%
bric_getAquisitionTime('/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/laser40_EM1000_exp100_001_stack.TIF', 'seconds')