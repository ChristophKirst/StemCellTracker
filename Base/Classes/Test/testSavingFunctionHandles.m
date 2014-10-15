% test classes

clear all 
close all

%%


t = Test


bkg = 0.1 * rand(3,4);
flt = rand(3,4);

%%

t.y = @(x)(correctFromBackgroudAndFlatField(x, bkg, flt))

%%

d = rand(3,4);
t.y(d)


%%

clear bkg
clear flt

%%

t.y

%%
save('test.mat', 't')


%%

clear all
close all

%%

load('test.mat')


t.y(rand(3,4))

%%
f = t.y

ff = functions(f)

ff.workspace{:}


