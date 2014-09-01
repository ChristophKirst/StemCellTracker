
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Qunatification Htt   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
initialize
bfinitialize

verbose = true;


addpath('./Scripts/User/Melissa')

%% Load Data

exp = Experiment('name', 'Melissa', 'description', 'some random segmentation for Melissa',...
                 'BaseDirectory', '/home/ckirst/Science/Projects/StemCells/Experiment/Data/Other/Melissa/', ...
                 'ImageDirectoryName', 'Results/',...
                 'ResultDirectoryName', 'Results/', ...
                 'ReadImageCommandFormat', 'load(''<file>'')',...
                 'ReadImageFileFormat', 'result_<sample>.mat');

exp.info()

%%
tags = exp.findImageTagRange

%% plot a smaple

filename = [exp.BaseDirectory, '/Images/140513_RUES2Ctrl_500um_11.zvi']
httplot(filename)


%% plot segmentation result

sample = 11
eval(exp.ReadImageCommand('sample', sample))





%%

filename = [exp.BaseDirectory, '/Images/140513_RUES2Q150_500um_32.zvi']
httplot(filename)




%% Counting on single samples -> histrogram


xvals = 400:40:16000;

control = [11];
%control = [11];

intC = {[],[],[]};
for sample = control
   eval(exp.ReadImageCommand('sample', sample))
   
   for c =1:3
      intC{c} = [intC{c} statsC{c}.MedianIntensity];
   end
end


figure(11); clf;
for c = 1:3
   subplot(3,1,c);
   {min(intC{c}), max(intC{c})}
   hist(intC{c}, xvals)
end


q150 = [31];

intQ = {[],[],[]};
for sample = q150
   eval(exp.ReadImageCommand('sample', sample))
   
   for c =1:3
      intQ{c} = [intQ{c} statsC{c}.MedianIntensity];
   end
end


figure(12); clf;
for c = 1:3
   subplot(3,1,c);
   {min(intQ{c}), max(intQ{c})}
   hist(intQ{c}, xvals)
end





%% pgenerate pdfs



xvals = 400:40:16000;

for sample = tags.sample

   control = [sample];
   %control = [11];
   
   intC = {[],[],[]};
   for sample = control
      eval(exp.ReadImageCommand('sample', sample))
      
      for c =1:3
         intC{c} = [intC{c} statsC{c}.MedianIntensity];
      end
   end
   
   
   h = figure(11); clf;
   for c = 1:3
      subplot(3,1,c);
      {min(intC{c}), max(intC{c})}
      hist(intC{c}, xvals)
   end
   
   
   print(h,'-dpdf', [exp.ResultDirectory, 'result_histogram_' num2str0(sample, 2) '.pdf']);

end








%% Control All


xvals = 400:40:16000;

control = [11 12 13 14 15 16 17 19 20 21];
%control = [11];


intC = {[],[],[]};
for sample = control
   eval(exp.ReadImageCommand('sample', sample))
   
   for c =1:3
      intC{c} = [intC{c} statsC{c}.MedianIntensity];
   end
end


h = figure(11); clf;
for c = 1:3
   subplot(3,1,c);
   {min(intC{c}), max(intC{c})}
   hist(intC{c}, xvals)
end


print(h,'-dpdf', [exp.ResultDirectory, 'result_histogram_ctrl.pdf']);




%% count 




%% Q150 All


q150 = [23 24 25 26 27 28 29 30 31 32];

intQ = {[],[],[]};
for sample = q150
   eval(exp.ReadImageCommand('sample', sample))
   
   for c =1:3
      intQ{c} = [intQ{c} statsC{c}.MedianIntensity];
   end
end


h = figure(12); clf;
for c = 1:3
   subplot(3,1,c);
   {min(intQ{c}), max(intQ{c})}
   hist(intQ{c}, xvals)
end

print(h,'-dpdf', [exp.ResultDirectory, 'result_histogram_q150.pdf']);


%% Total numbers


ctotal = length(intC{1})
cq150 = length(intQ{1})



%% thresholding

clc

thbra = 1500;
thsox2 = 1800;
thcdx2 = 3000;

ids = (intC{1} > thcdx2);
ccdx2 = total(ids)


ids = (intC{2} > thbra);
cbra = total(ids)

ids = (intC{3} > thsox2);
csox2 = total(ids)



ids = (intQ{1} > thcdx2);
qcdx2 = total(ids)

ids = (intQ{2} > thbra);
qbra = total(ids)

ids = (intQ{3} > thsox2);
qsox2 = total(ids)



