%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test File Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all


%% dirr

dd = dirr('*.tif')

dd = dirr('*')

dd = dirr('*', 'directories', true, 'files', false)


%%

dd = dirr('*/Test', 'directories', true)


%%

dd = dirr('*.tif', 'recursive', true)





%% Infer Tagformat from file list




%% FileHandler



fh = FileHandler('BaseDirectory',          './Test/Data/Experiment', ...
                 'ImageDirectoryName',     '../../Images/hESCells_Folder',...
                 'ResultDirectoryName',    '.', ...
                 'ReadImageCommandFormat', 'imread(''<file>'')',...
                 'ReadImageFileFormat',    'f<field>_t<time>/f<field>_t<time,2>_<str,1,s><zpos,2>.tif');

fh.info()


%% Get file tags


tags = fh.findImageTags

fh.findImageTagRange



%% Infer Tagformat form file list


fns = fh.findImageFiles;

tagformat(fns)




%% Test some Regexps

fname = 'f1_t02.tif';
fname2 = 'f1_t01.tif';
expr =  tagformat2regexp('f<f>_t<t,2>.tif')


fname = 'f1_t1/f1_t02.tif';
fname2 = 'f1_t1/f1_t01.tif';
expr =  tagformat2regexp('f<f>_t<t>/f<f>_t<t,2>.tif')


regexp(fname, expr)
regexp(fname2, expr)
