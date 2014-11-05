%%%%%%%%%%%%%%%%%%%%
%%% Test Python %%%%
%%%%%%%%%%%%%%%%%%%%

initialize

clc
clear all

%%


m = 1:6;
m = reshape(m, 2,3)
p = numpyFromMat(m)
m

%%

clc
m2 = numpyToMat(p);

m - m2


%% large data sets

clear all
clc


%%

m = rand(2000, 2000, 15);


%%
tic
p = numpyFromMat(m)
toc


%%

p.shape

%%

tic
m2 = numpyToMat(p);
toc
