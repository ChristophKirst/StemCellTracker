

%% Load Data


clear all
clear java

initialize
bfinitialize

filename = '/home/ckirst/Science/Projects/StemCells/Experiment/Other/Wnt/wnt_clone8_again_feb09.lif';
seriesid = 1;
datalif = imread_bf(filename, struct('series', seriesid, 'channel', []));

%%


for t = 1:5
 
   stack = datalif(:,:,:,1,t);
   sz = size(stack,3);
   
   for s = 1:sz
      imwrite(stack(:,:,s), ['/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/mESCells_Wnt/wnt_t' num2str0(t,2) '_z' num2str0(s,2), '.tif']);
   end

end
   


