%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test OverlapRegoins  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


isizes = repmat({[10,10]},1,3)
shifts = {[0,0], [5,0], [3,5]}; %, [0,5], [3,3]}

clc
[reg, ids] = stitchImagesOverlapRegions(shifts, isizes);


%%
var2char(reg)


img = zeros(15);

for i = 1:length(reg)
   te = reg{i};
   img(te(1):te(2), (te(3)):(te(4))) = i;
end

figure(1); clf
implot(img);


unique(img(:))

   