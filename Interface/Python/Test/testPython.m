%%%%%%%%%%%%%%%%%%%%
%%% Test Python %%%%
%%%%%%%%%%%%%%%%%%%%

initialize

%%


m = rand(5,6);
p = numpyFromMat(m)

%%

clc
m2 = numpyToMat(p);

m - m2


%% large data sets

m = rand(10000);


%%
tic
p = numpyFromMat(m)
toc

%%

tic
m2 = numpyToMat(p);
toc
