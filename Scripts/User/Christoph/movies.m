%%%%%%%%%%%%%%%%%%%%%
%%% Create Movies %%%
%%%%%%%%%%%%%%%%%%%%%



%% long term imaging
times = {'14_56', '15_16', '15_43', '16_09', '16_35'};
stage = 's1';
rpath = '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/';

fns = {};
for i = 1:length(times)
   fn = dir([rpath times{i} '/' stage '/*.TIF']);    
   fns = [fns, cellfun(@(x) fullfile(rpath, times{i}, stage, x), {fn.name}, 'UniformOutput', false)];
end

length(fns)


img = imread(fns{7});
max(img(:))

%%

mov = images2movie(fns,...
           ['/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/movie_' stage '.avi'],... 
           'norm', 2400);



        
        
%% search for large scale events

rpath = '/var/run/media/ckirst/ChristophsData/Science/Projects/StemCells/Experiment/Calcium/hES_GCamp6/2014_04_04/12_43_flash_search/';

fns = dir(fullfile(rpath, 's1', '*.TIF'));
fns = cellfun(@(x) fullfile(rpath, 's1', x), {fns.name}, 'UniformOutput', false);

length(fns)

times = mmaquisitiontime(fns);
